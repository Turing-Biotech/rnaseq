#!/usr/bin/env python3
"""Render a Chinese (ZH) fusion QC report for nf-core/rnaseq.

Reads raw tool outputs (not MultiQC HTML) from a finished pipeline output
directory, plus optional methylation contract files, and renders one
self-contained HTML report (plus best-effort PDF via wkhtmltopdf).

Design notes:
  * Every section degrades to an explicit "NA / 暂无数据" block if its source
    files are missing - the report always renders.
  * Methylation sections only appear with data when --meth-dir contains the
    vendor contract files (FilterSummary.xls etc., which are actually TSV).
  * No hard dependency on pandas: stdlib csv/json/zipfile only (jinja2 + pyyaml).
"""

import argparse
import csv
import glob
import html as htmllib
import os
import re
import shutil
import subprocess
import sys
import zipfile

try:
    import yaml
except ImportError:  # pragma: no cover
    yaml = None

try:
    import jinja2
except ImportError:  # pragma: no cover
    jinja2 = None

NA = "NA"

# --------------------------------------------------------------------------
# small helpers
# --------------------------------------------------------------------------

def esc(v):
    return htmllib.escape(str(v), quote=False)


def to_num(s):
    """Loose numeric parsing (handles '12,345', '98.0', '1.2e5', 'NA')."""
    if s is None:
        return None
    s = str(s).strip().replace(",", "").rstrip("%")
    m = re.search(r"-?\d+\.?\d*(?:[eE][-+]?\d+)?", s)
    if not m:
        return None
    try:
        return float(m.group(0))
    except ValueError:
        return None


def fmt_int(v):
    if v is None:
        return NA
    return "{:,.0f}".format(v)


def fmt_num(v, nd=2):
    if v is None:
        return NA
    if abs(v) >= 100:
        return "{:.1f}".format(v)
    return "{:.{n}f}".format(v, n=nd)


def fmt_bp(v):
    if v is None:
        return NA
    for unit, lim in (("Gb", 1e9), ("Mb", 1e6)):
        if v >= lim:
            return "{:.2f} {}".format(v / lim, unit)
    return "{:,.0f} bp".format(v)


def read_tsv(path):
    """Read a .txt/.xls-that-is-TSV with no strict header assumption.

    Returns (header, rows) with rows as list[str].
    """
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        text = fh.read()
    if not text.strip():
        return [], []
    lines = [ln.rstrip("\n") for ln in text.splitlines() if ln.strip()]
    delim = "\t" if "\t" in lines[0] else ("," if "," in lines[0] else "\t")
    rdr = csv.reader(lines, delimiter=delim)
    rows = [r for r in rdr if any(c.strip() for c in r)]
    if not rows:
        return [], []
    return rows[0], rows[1:]


def html_table(header, rows, na_fill=True):
    out = ['<table>']
    if header:
        out.append("<tr>" + "".join("<th>{}</th>".format(esc(h)) for h in header) + "</tr>")
    for r in rows:
        cells = []
        for i, c in enumerate(r):
            if c is None or str(c).strip() == "":
                c = NA if (na_fill and i > 0) else ""
            cells.append("<td>{}</td>".format(esc(c)))
        if cells:
            out.append("<tr>" + "".join(cells) + "</tr>")
    out.append("</table>")
    return "".join(out)


def na_block(reason="该章节暂无数据（上游流程未产出对应结果，或参数未提供）"):
    return "<p class='na'>{}</p>".format(esc(reason))


def find_files(root, *suffix_patterns):
    """Glob files under root whose name matches any regex (case-insensitive)."""
    if not root or not os.path.isdir(root):
        return []
    hits = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            if fn.startswith("~$") or fn.startswith("."):
                continue
            for pat in suffix_patterns:
                if re.search(pat, fn, re.IGNORECASE):
                    p = os.path.join(dirpath, fn)
                    hits.append(p)
                    break
    return sorted(hits)

# --------------------------------------------------------------------------
# samplesheet
# --------------------------------------------------------------------------

def load_samplesheet(path):
    """Return (groups, individuals, raw_header, raw_rows).

    groups: OrderedDict group -> [individual labels]
    individuals: list of (label, group, fastq1, strandedness)
    """
    groups = {}
    individuals = []
    if not path or not os.path.isfile(path):
        return groups, individuals, None, []
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        rdr = csv.DictReader(fh)
        header = rdr.fieldnames or []
        for row in rdr:
            sample = (row.get("sample") or "").strip()
            if not sample:
                continue
            gcol = row.get("group")
            group = gcol.strip() if gcol and gcol.strip() else sample
            fq1 = row.get("fastq_1") or ""
            label = os.path.basename(fq1)
            label = re.sub(r"\.(fastq|fq)\.gz$", "", label)
            label = re.sub(r"[-_.]?(R|r)[12](_\d+)?(.fastq|.fq)?\.gz$", "", label)
            label = label or sample
            individuals.append((label, group, fq1, (row.get("strandedness") or "").strip()))
            groups.setdefault(group, [])
            if label not in groups[group]:
                groups[group].append(label)
    return groups, individuals, header, list(header)


# --------------------------------------------------------------------------
# RNA-seq section loaders
# --------------------------------------------------------------------------

def section_groups(groups, meth_group_plan_rows=None):
    """s1_2_1: 实验分组"""
    if not groups:
        return na_block("samplesheet 中未解析到样本分组")
    header = ["组别", "样本数", "样本列表"]
    rows = [[g, len(v), "、".join(v)] for g, v in groups.items()]
    note = ""
    if meth_group_plan_rows:
        gp_header, gp_rows = meth_group_plan_rows[0], meth_group_plan_rows[1]
        if gp_rows:
            header, rows = gp_header, gp_rows
            note = "<p class='na'>上表来自甲基化端 GroupPlan（含全部个体）；RNA-seq 端按上表组别聚合。</p>"
    return html_table(header, rows) + note


def section_deg_comparison(deg_comparison, groups, meth_dmr_rows=None):
    """s1_2_2: 差异分析对比"""
    blocks = []
    rows = []
    if deg_comparison:
        rows.append([deg_comparison, "DESeq2 (DESEQ2_QC)"])
    elif len(groups) == 2:
        a, b = list(groups.keys())
        rows.append(["{} vs {}".format(b, a), "DESeq2 (DESEQ2_QC)"])
    if rows:
        blocks.append(html_table(["差异比较", "分析工具"], rows))
    if meth_dmr_rows and meth_dmr_rows[1]:
        dmr_header, dmr_body = meth_dmr_rows[0], meth_dmr_rows[1]
        blocks.append("<p class='cap'>甲基化端 DMR 对比 (DMRPlan, 分表)</p>"
                      + html_table(dmr_header, dmr_body))
    if not blocks:
        return na_block("未提供 --deg-comparison，且样本不足以推断对比设计")
    return "".join(blocks)


def section_reference(genome, gtf, project_id, meta_extra):
    """s2: 参考信息"""
    rows = []
    if genome:
        rows.append(["参考基因组", genome])
    if gtf:
        rows.append(["GTF 注释", gtf])
    if project_id:
        rows.append(["项目编号", project_id])
    for k, v in (meta_extra or {}).items():
        if v:
            rows.append([k, v])
    if not rows:
        return na_block()
    return html_table([], rows)


def _parse_trimgalore(path):
    """Parse a Trim Galore report (2.1.x cutadapt-style or classic)."""
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    except OSError:
        return None
    out = {}
    m = re.search(r"Total reads processed:\s*([\d,]+)", text)
    if m:
        out["raw_reads"] = to_num(m.group(1))
    m = re.search(r"Reads written \(passing filters\):\s*([\d,]+)\s*\(([\d.]+)%\)", text)
    if m:
        out["clean_reads"], out["clean_pct"] = to_num(m.group(1)), to_num(m.group(2))
    m = re.search(r"Total basepairs processed:\s*([\d,]+)", text)
    if m:
        out["raw_bp"] = to_num(m.group(1))
    m = re.search(r"Total written \(filtered\):\s*([\d,]+)\s*bp\s*\(([\d.]+)%\)", text)
    if m:
        out["clean_bp"], out["clean_bp_pct"] = to_num(m.group(1)), to_num(m.group(2))
    # classic trimgalore format fallback
    if "raw_reads" not in out:
        m = re.search(r"Input\s*:\s*([\d,]+)\s*reads", text)
        if m:
            out["raw_reads"] = to_num(m.group(1))
        m = re.search(r"Clean\s*:\s*([\d,]+)\s*reads.*?\(([\d.]+)%\)", text, re.S)
        if m:
            out["clean_reads"], out["clean_pct"] = to_num(m.group(1)), to_num(m.group(2))
    for q, key in (("Q20", "q20"), ("Q30", "q30")):
        m = re.search(r"{0}\s*:\s*([\d.]+)%".format(q), text)
        if m:
            out[key] = to_num(m.group(1))
    return out or None


def _parse_fastp_json(path):
    try:
        import json as _json
        with open(path, "r", encoding="utf-8") as fh:
            d = _json.load(fh)
    except Exception:
        return None
    out = {}
    for r in ("before_filtering", "after_filtering"):
        b = d.get(r) or {}
        out["raw_reads"] = to_num(b.get("total_reads")) if r == "before_filtering" else out.get("raw_reads")
        out["clean_reads"] = to_num(b.get("filtered_reads")) if r == "after_filtering" else out.get("clean_reads")
        out["raw_bp"] = to_num(b.get("total_bases")) if r == "before_filtering" else out.get("raw_bp")
        out["clean_bp"] = to_num(b.get("filtered_bases")) if r == "after_filtering" else out.get("clean_bp")
    pf = d.get("filtering_result") or {}
    out["clean_pct"] = (100.0 * pf.get("filtered_reads", 0) / pf.get("total_reads", 1)) if pf.get("total_reads") else None
    mtr = d.get("basic_raw") or {}
    out["q20"] = to_num(mtr.get("q20_rate"))
    mtr2 = d.get("basic_rejected") or d.get("filtering_result") or {}
    if "q20" not in out:
        q20 = d.get("filtering_result") or {}
        out["q20"] = to_num(q20.get("q20_rate"))
    return out or None


def section_trimming(outdir, groups):
    """s3_1: 过滤前后对比 (RNA 端; 甲基化端由 loader 追加)"""
    stats = {}
    for g in groups:
        # trimgalore
        cands = glob.glob(os.path.join(outdir, "trimgalore", g + "_*trimming_report.txt"))
        r1 = [p for p in cands if re.search(r"[-_.]1[-_.](fastq|fq)\.gz_trimming_report\.txt$", p)]
        cand = (r1 or cands)
        data = None
        if cand:
            data = _parse_trimgalore(cand[0])
        if data is None:
            fjs = glob.glob(os.path.join(outdir, "fastp", "*.json"))
            fjs = [f for f in fjs if g in os.path.basename(f)]
            if fjs:
                data = _parse_fastp_json(fjs[0])
        stats[g] = data or {}
    if not any(stats.values()):
        return na_block("未找到 trimgalore/fastp 过滤报告")
    header = ["组别", "下机量 (reads)", "清洁度 (reads)", "清洁度占比 (%)", "下机量 (碱基)", "清洁度 (碱基)", "Q20 (%)", "Q30 (%)"]
    rows = []
    for g in groups:
        d = stats[g]
        rows.append([
            g,
            fmt_int(d.get("raw_reads")),
            fmt_int(d.get("clean_reads")),
            fmt_num(d.get("clean_pct")),
            fmt_bp(d.get("raw_bp")),
            fmt_bp(d.get("clean_bp")),
            fmt_num(d.get("q20")),
            fmt_num(d.get("q30")),
        ])
    return html_table(header, rows)


def _fastqc_zip_stats(zip_path):
    """Parse fastqc_data.txt (FastQC 0.12 tab format, with legacy fallbacks)."""
    try:
        with zipfile.ZipFile(zip_path) as zf:
            names = [n for n in zf.namelist() if n.endswith("fastqc_data.txt")]
            if not names:
                return {}
            text = zf.read(names[0]).decode("utf-8", errors="replace")
    except Exception:
        return {}
    out = {}
    m = re.search(r"Total Sequences\t([\d,\.]+)", text)
    if m:
        out["nseq"] = to_num(m.group(1))
    m = re.search(r"%GC\t?([\d\.]+)", text)
    if m:
        out["gc"] = to_num(m.group(1))
    m = re.search(r"#Total Deduplicated Percentage\t?([\d\.]+)", text)
    if m:
        v = to_num(m.group(1))
        if v is not None:
            out["dup_pct"] = 100.0 - v
    sec = re.search(r">>Per sequence quality scores.*?\n(.*?)>>END_MODULE", text, re.S)
    if sec:
        total = 0.0
        wsum = 0.0
        q30c = 0.0
        for ln in sec.group(1).splitlines():
            ln = ln.strip()
            if not ln or ln.startswith("#"):
                continue
            parts = ln.split("\t")
            if len(parts) < 2:
                continue
            v = to_num(parts[0])
            c = to_num(parts[1])
            if v is None or c is None:
                continue
            total += c
            wsum += v * c
            if int(v) == 30:
                q30c = c
        if total:
            out["avg_seq_q"] = wsum / total
            out["q30_pct"] = 100.0 * q30c / total
    # legacy (space+colon) fallbacks
    m = re.search(r"Average Quality Score:\s*([\d\.]+)", text)
    if m and "avg_seq_q" not in out:
        out["avg_seq_q"] = to_num(m.group(1))
    m = re.search(r"%GC:\s*([\d\.]+)", text)
    if m and "gc" not in out:
        out["gc"] = to_num(m.group(1))
    return out


def section_fastqc(outdir, groups):
    """s3_2: FastQC 质量 (raw 阶段, 每组取 R1)"""
    rows = []
    header = ["组别", "reads 数", "平均质量分 (per-read)", "Q30 占比 (%)", "%GC", "重复率 (%)"]
    for g in groups:
        zips = glob.glob(os.path.join(outdir, "fastqc", "raw", "**", "*.zip"), recursive=True)
        zips = [z for z in zips if re.search(
            re.escape(g) + r"[-_](raw[-_])?1(-val-1)?_fastqc\.zip$", os.path.basename(z))]
        if not zips:
            continue
        d = _fastqc_zip_stats(zips[0])
        if not d:
            continue
        rows.append([
            g,
            fmt_int(d.get("nseq")),
            fmt_num(d.get("avg_seq_q"), 1),
            fmt_num(d.get("q30_pct"), 1),
            fmt_num(d.get("gc"), 1),
            fmt_num(d.get("dup_pct"), 1),
        ])
    if not rows:
        return na_block("未找到 FastQC 报告 (fastqc zip)")
    return html_table(header, rows)


def _parse_star_log(path):
    try:
        with open(path, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    except OSError:
        return None
    out = {}
    pats = {
        "total": r"Number of input reads\s*\|\s*(.+)",
        "uniq_pct": r"Uniquely mapped reads %\s*\|\s*([\d.]+)%",
        "multi_pct": r"% of reads mapped to multiple loci\s*\|\s*([\d.]+)%",
        "dups_pct": r"Duplicated reads %\s*\|\s*([\d.]+)%",
        "conv_pct": r"Conversion rate\s*\|\s*([\d.]+)%",
    }
    for k, p in pats.items():
        m = re.search(p, text)
        if m:
            out[k] = to_num(m.group(1))
    # combined unmapped
    unmapped = None
    for key in ("too short", "other"):
        m = re.search(r"% of reads unmapped: " + key + r"\s*\|\s*([\d.]+)%", text)
        if m:
            unmapped = (unmapped or 0) + to_num(m.group(1))
    if unmapped is not None:
        out["unmapped_pct"] = unmapped
    return out or None


def section_alignment(outdir, groups):
    """s4_1: 比对统计 (STAR)"""
    rows = []
    for g in groups:
        logs = glob.glob(os.path.join(outdir, "star_salmon", "log", g + ".Log.final.out")) + \
               glob.glob(os.path.join(outdir, "star", "log", g + ".Log.final.out"))
        if logs:
            d = _parse_star_log(logs[0])
            if d:
                rows.append([
                    g,
                    fmt_int(d.get("total")),
                    fmt_num(d.get("uniq_pct")),
                    fmt_num(d.get("multi_pct")),
                    fmt_num(d.get("unmapped_pct")),
                    fmt_num(d.get("dups_pct")),
                    fmt_num(d.get("conv_pct")),
                ])
    if not rows:
        return na_block("未找到 STAR Log.final.out 日志")
    header = ["组别", "输入 reads", "唯一比对 (%)", "多比对 (%)", "未比对 (%)", "重复率 (%)", "转换率 (%)"]
    return html_table(header, rows)


def _parse_salmon_quant(gdir):
    out = {}
    gsf = os.path.join(gdir, "quant.genes.sf")
    tsf = os.path.join(gdir, "quant.sf")
    if os.path.isfile(gsf):
        header, rows = read_tsv(gsf)
        hi = {h.strip().lower(): i for i, h in enumerate(header)}
        tpi = hi.get("tpm")
        ci = hi.get("estimated_count", hi.get("count"))
        n = 0
        if tpi is not None:
            for r in rows:
                if ci is not None and ci < len(r) and (to_num(r[ci]) or 0) > 0:
                    n += 1
                elif tpi < len(r) and (to_num(r[tpi]) or 0) > 0:
                    n += 1
        out["detected_genes"] = n
    elif os.path.isfile(tsf):
        header, rows = read_tsv(tsf)
        hi = {h.strip().lower(): i for i, h in enumerate(header)}
        tpi = hi.get("tpm")
        n = 0
        tot = 0.0
        for r in rows:
            v = to_num(r[tpi]) if tpi is not None and tpi < len(r) else None
            if v:
                tot += v
                if v > 0.5:
                    n += 1
        out["detected_transcripts"] = n
        out["total_tpm"] = tot
    return out


def section_quantification(outdir, groups):
    """s6: 转录组定量"""
    rows = []
    for g in groups:
        gdir = os.path.join(outdir, "star_salmon", g)
        if not os.path.isdir(gdir):
            gdir = os.path.join(outdir, "salmon", g)
        d = _parse_salmon_quant(gdir) if os.path.isdir(gdir) else {}
        if not d:
            continue
        rows.append([g, fmt_int(d.get("detected_genes") or d.get("detected_transcripts"))])
    if not rows:
        return na_block("未找到定量结果 (salmon quant.sf / quant.genes.sf)")
    header = ["组别", "检出基因数 (count > 0)"]
    html = html_table(header, rows)
    merged = glob.glob(os.path.join(outdir, "star_salmon", "*merged*"))
    if merged:
        html += "<p class='na'>定量矩阵文件: {}</p>".format(esc(", ".join(os.path.relpath(m, outdir) for m in merged[:6])))
    return html


# --------------------------------------------------------------------------
# DeSeq2 QC
# --------------------------------------------------------------------------

def section_deseq2(deseq2_dir):
    """s7_1: PCA 与样本间距离 (from DESEQ2_QC mqc TSVs)"""
    blocks = []
    if deseq2_dir and os.path.isdir(deseq2_dir):
        for pat, cap in (("*pca.vals_mqc.tsv", "PCA 主坐标"), ("*dists_mqc.tsv", "样本间距离矩阵")):
            for f in sorted(glob.glob(os.path.join(deseq2_dir, pat))):
                header, rows = read_tsv(f)
                if rows:
                    blocks.append("<p class='cap'>" + esc(cap) + "</p>" + html_table(header, rows))
                else:
                    blocks.append("<p class='cap'>" + esc(cap) + "</p><p class='na'>文件为空</p>")
    if not blocks:
        return na_block("未找到 DeSeq2 QC 表格 (PCA / 距离矩阵)")
    return "".join(blocks)


# --------------------------------------------------------------------------
# methylation-side loaders (vendor contract files, actually TSV)
# --------------------------------------------------------------------------

def _meth_pick(meth_files, *patterns):
    for pat in patterns:
        for p in meth_files:
            if re.search(pat, p, re.IGNORECASE):
                return p
    return None


def section_meth_filter(meth_files, groups):
    """甲基化端过滤统计 (追加到 3.1); 返回 (table_html or None, sample_names)"""
    f = _meth_pick(meth_files, r"FilterSummary")
    if not f:
        return None, []
    header, rows = read_tsv(f)
    if not rows:
        return None, []
    names = []
    for r in rows:
        if r and r[0].strip():
            names.append(r[0].strip())
    return html_table(header, rows), names


def section_meth_alignment(meth_files):
    f = _meth_pick(meth_files, r"GenomeMapping|GenomeMappingSummary")
    if not f:
        return na_block("甲基化端未提供 GenomeMappingSummary")
    header, rows = read_tsv(f)
    if not rows:
        return na_block("GenomeMappingSummary 为空")
    return html_table(header, rows)


def section_csites(meth_files):
    """s4_2: mC 覆盖 (pivot: rows=samples, cols=coverage bins)"""
    per = {}
    bins_order = []
    for p in meth_files:
        if not re.search(r"CsitesCoverage|Csites_coverage|csites", p, re.IGNORECASE):
            continue
        header, rows = read_tsv(p)
        if not rows:
            continue
        sample = os.path.basename(p).rsplit(".CsitesCoverage", 1)[0]
        sample = re.sub(r"\.(xls|tsv|txt|csv)$", "", sample)
        d = {}
        for r in rows:
            if len(r) >= 2:
                d[r[0].strip()] = r[1].strip()
                if r[0].strip() not in bins_order:
                    bins_order.append(r[0].strip())
        per[sample] = d
    if not per:
        return na_block("未找到 CsitesCoverage 覆盖统计文件")
    header = ["样本"] + bins_order
    rows = [[s] + [d.get(b, "") for b in bins_order] for s, d in per.items()]
    return "<p class='cap'>mC 覆盖度分布（百分比，例示）</p>" + html_table(header, rows)


def _mesh_stat_table(meth_files, *patterns):
    f = _meth_pick(meth_files, *patterns)
    if not f:
        return na_block("未找到对应统计文件")
    header, rows = read_tsv(f)
    if not rows:
        return na_block("统计文件为空")
    return html_table(header, rows)


def _mesh_stat_tables(meth_files, pattern, cap):
    """Render one table per matching stat file (vendor files are per-sample)."""
    hits = [p for p in meth_files if re.search(pattern, p, re.IGNORECASE)]
    if not hits:
        return na_block("未找到对应统计文件")
    blocks = []
    for f in hits:
        header, rows = read_tsv(f)
        if not rows:
            continue
        sample = os.path.basename(f).rsplit(".chr_methy_stat", 1)[0]
        sample = re.sub(r"\.(element_methy_stat|chr_methy_stat|number_stat)\.xls$", "", sample)
        sample = sample.replace("allsample", "全部样本")
        blocks.append("<p class='cap'>{0} · {1}</p>".format(esc(cap), esc(sample)) + html_table(header, rows))
    if not blocks:
        return na_block("统计文件为空")
    return "".join(blocks)


def section_meth_chr(meth_files):
    """s5_1: 染色体水平分布"""
    return _mesh_stat_tables(meth_files, r"chr_methy_stat|chr\.?number_stat", "染色体水平")


def section_meth_element(meth_files):
    """s5_2: 元件水平分布"""
    return _mesh_stat_tables(meth_files, r"element_methy_stat|element\.?number_stat", "元件水平")


def section_meth_number(meth_files):
    """s5_3: 不同甲基化水平分类统计"""
    return _mesh_stat_tables(meth_files, r"number_stat", "甲基化片段水平")


def section_meth_presence(meth_files):
    if meth_files:
        return "<p class='na'>检测到 {} 个甲基化端契约文件，相应章节已填充数据。</p>".format(len(meth_files))
    return na_block("未提供甲基化端结果目录（--meth-result-dir）；本运行仅包含 RNA-seq 报告内容。")


# --------------------------------------------------------------------------
# versions + methods
# --------------------------------------------------------------------------

def section_versions(versions_path):
    if not versions_path or not os.path.isfile(versions_path) or yaml is None:
        return na_block("未找到软件版本清单 (versions.yml)")
    try:
        with open(versions_path, "r", encoding="utf-8") as fh:
            data = yaml.safe_load(fh) or {}
    except Exception as exc:
        return na_block("versions.yml 解析失败: " + str(exc))
    rows = []
    for proc, mods in sorted(data.items()):
        if not isinstance(mods, dict):
            continue
        for mod, ver in sorted(mods.items()):
            rows.append([proc.replace("_", " "), mod, ver])
    if not rows:
        return na_block("versions.yml 中无条目")
    return html_table(["流程步骤", "组件", "版本"], rows)


def section_methods(command_line, genome, title, company):
    text = (
        "本报告由 nf-core/rnaseq 流程自动运行生成。Raw data 经 adapter 修剪与质量过滤后，"
        "使用 STAR 将 reads 比对到参考基因组" + (" ({})".format(esc(genome)) if genome else "") +
        "，比对后的 reads 使用 salmon 进行转录本定量，定量矩阵使用 R/DESeq2 进行 PCA 等 QC。"
        + ("甲基化端结果（Read 过滤、比对与甲基化水平统计）来自配套甲基化分析流程的标准化输出。" if True else "")
        + " 完整参数与命令见下。报告生成: " + esc(company or "nf-core/rnaseq") + "。"
    )
    html = "<p>{}</p>".format(text)
    if command_line:
        html += "<pre><code>" + esc(command_line) + "</code></pre>"
    return html


def section_degenerics(deg_table):
    """s7_2: 差异表达概况 (external DEG csv, optional)"""
    if not deg_table or not os.path.isfile(deg_table):
        return na_block("未提供 --deg-table，本小节默认跳过")
    try:
        header, rows = read_tsv(deg_table)
    except Exception as exc:
        return na_block("DEG 表读取失败: " + str(exc))
    if not rows:
        return na_block("DEG 表为空")
    # count up/down
    hi = {h.strip().lower(): i for i, h in enumerate(header)}
    lfc_i = None
    pval_i = None
    for k in ("log2fc", "log2 fold change", "logfc", "log2foldchange"):
        if k in hi:
            lfc_i = hi[k]
            break
    for k in ("padj", "pvalue", "p value", "adjpvalue"):
        if k in hi:
            pval_i = hi[k]
            break
    up = down = 0
    for r in rows:
        try:
            lfc = to_num(r[lfc_i]) if lfc_i is not None and lfc_i < len(r) else None
            pv = to_num(r[pval_i]) if pval_i is not None and pval_i < len(r) else None
            if lfc is not None and (pv is None or pv < 0.05):
                if lfc > 0:
                    up += 1
                elif lfc < 0:
                    down += 1
        except Exception:
            pass
    note = "<p class='cap'>共检出差异表达: 上调 {} 个, 下调 {} 个 (|padj| &lt; 0.05)。下表展示前 30 行 (共 {} 行)。</p>".format(up, down, len(rows))
    return note + html_table(header, rows[:30])


def section_abstract(groups, individuals, total_stats):
    parts = []
    parts.append("本项目为转录组 RNA-seq 实验分析。")
    if groups:
        parts.append("共 {} 个组别 ({}), {} 个样本。".format(
            len(groups), "、".join(groups.keys()), len(individuals)))
    if total_stats:
        parts.append("下机总量约 {} reads, 过滤后清洁度 {}。".format(
            fmt_int(total_stats.get("total_reads")),
            fmt_num(total_stats.get("overall_clean_pct")) + "%" if total_stats.get("overall_clean_pct") is not None else NA))
    parts.append("流程输出包括 Read 过滤、比对、定量与表达量 QC; 各节详见下文。")
    return "<p>{}</p>".format("".join(parts))


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------

def load_terms(path):
    """Section registry from report_zh_terms.yml; falls back to built-in."""
    default = {
        "s1": (1, "项目概述"), "s1_1": (2, "摘要"), "s1_2": (2, "实验设计"),
        "s1_2_1": (3, "实验分组"), "s1_2_2": (3, "差异分析对比"),
        "s2": (1, "参考信息"), "s3": (1, "Read 过滤"), "s3_1": (2, "过滤前后对比"),
        "s3_2": (2, "FastQC 质量"), "s4": (1, "Read 比对"), "s4_1": (2, "比对统计"),
        "s4_2": (2, "mC 覆盖"), "s5": (1, "甲基化水平统计"), "s5_1": (2, "染色体水平分布"),
        "s5_2": (2, "元件水平分布"), "s5_3": (2, "不同甲基化水平分类统计"),
        "s6": (1, "转录组定量"), "s7": (1, "表达量 QC"), "s7_1": (2, "PCA 与样本间距离"),
        "s7_2": (2, "差异表达概况"), "s8": (1, "软件版本与方法"), "s8_1": (2, "软件版本"),
        "s8_2": (2, "分析方法"),
    }
    if path and os.path.isfile(path) and yaml is not None:
        try:
            with open(path, "r", encoding="utf-8") as fh:
                d = yaml.safe_load(fh) or {}
            secs = d.get("sections") or {}
            for k, v in secs.items():
                if isinstance(v, dict) and "level" in v and "title" in v:
                    default[k] = (int(v["level"]), str(v["title"]))
        except Exception:
            pass
    return default


def make_pdf(html_path, pdf_path):
    exe = shutil.which("wkhtmltopdf")
    if not exe:
        return False, "wkhtmltopdf not available"
    cmd = [exe, "--quiet", "--encoding", "utf-8", "--enable-local-file-access",
           "--page-size", "A4", "--margin-top", "18mm", "--margin-bottom", "18mm",
           "--margin-left", "16mm", "--margin-right", "16mm", html_path, pdf_path]
    try:
        r = subprocess.run(cmd, timeout=900, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if r.returncode == 0 and os.path.isfile(pdf_path):
            return True, "ok"
        return False, "wkhtmltopdf exit {}: {}".format(r.returncode, (r.stderr or b"")[:300].decode("utf-8", "replace"))
    except Exception as exc:
        return False, str(exc)


def _unquote(v):
    v = (v or "").strip()
    if v.startswith("'") and v.endswith("'"):
        return v[1:-1].replace("'\\''", "'")
    return v


CSS = """
body { font-family: "Noto Sans CJK SC","Noto Sans CJK","WenQuanYi Zen Hei","Microsoft YaHei",sans-serif;
       font-size: 10pt; color: #222; line-height: 1.5; }
h1 { font-size: 15pt; color: #1f3864; border-bottom: 2px solid #2f5597; padding-bottom: 3px; margin-top: 22px; }
h2 { font-size: 12pt; color: #2f5597; margin-top: 16px; }
h3 { font-size: 10.5pt; color: #44546a; margin-top: 12px; }
table { border-collapse: collapse; width: 100%; font-size: 8pt; margin: 8px 0; }
th { background: #2f5597; color: #fff; padding: 3px 5px; border: 1px solid #9db2d5; font-weight: 600; }
td { border: 1px solid #b8c4d9; padding: 2.5px 5px; text-align: center; word-break: break-all; }
tr:nth-child(even) td { background: #f2f5fa; }
.na { color: #8a8a8a; font-size: 9pt; }
.cap { color: #44546a; font-size: 9pt; font-weight: 600; margin: 8px 0 2px; }
pre { background: #f5f6f8; border: 1px solid #d7dbe2; padding: 8px; font-size: 8pt; white-space: pre-wrap; word-break: break-all; }
.cover { text-align: center; margin-top: 110px; page-break-after: always; }
.cover img { max-width: 320px; margin-bottom: 26px; }
.cover .company { font-size: 15pt; color: #1f3864; letter-spacing: 3px; margin-bottom: 60px; }
.cover .title { font-size: 26pt; font-weight: 700; color: #1f3864; margin-bottom: 90px; }
.cover .dated { font-size: 11pt; color: #555; }
.toc { page-break-after: always; }
.toc h1 { border: none; }
.toc ul { list-style: none; padding-left: 0; }
.toc li { margin: 3px 0; }
.toc a { text-decoration: none; color: #222; }
.l2 { padding-left: 1.6em; color: #44546a; }
.l3 { padding-left: 3.2em; color: #66707e; font-size: 9pt; }
@page { @bottom-center { content: counter(page) " / " counter(pages); font-size: 8pt; color: #888; } }
"""

TEMPLATE = r"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>{css}</style>
</head>
<body>
<div class="cover">
{logo_html}
<div class="company">{company}</div>
<div class="title">{title}</div>
<div class="dated">{date}</div>
</div>
<div class="toc">
<h1>目录</h1>
<ul>
{toc}
</ul>
</div>
{sections}
</body>
</html>
"""


def render_report(args):
    sections_cfg = load_terms(args.terms)
    outdir = args.results
    groups, individuals, _h, _r = load_samplesheet(args.samplesheet)

    # aggregate trim stats for the abstract
    total_stats = {}
    if outdir and os.path.isdir(outdir):
        raw = clean = 0
        for g in groups:
            for p in glob.glob(os.path.join(outdir, "trimgalore", g + "_*_trimming_report.txt")):
                if re.search(r"[-_.]1[-_.](fastq|fq)\.gz_trimming_report\.txt$", p):
                    d = _parse_trimgalore(p) or {}
                    raw += d.get("raw_reads") or 0
                    clean += d.get("clean_reads") or 0
                    break
        if raw:
            total_stats = {"total_reads": raw, "overall_clean_pct": 100.0 * clean / raw}

    meth_files = args.meth_dir and [f for f in (os.listdir(args.meth_dir) if os.path.isdir(args.meth_dir) else [])] or []
    meth_files = [os.path.join(args.meth_dir, f) for f in meth_files] if meth_files else []

    # meth group/dmr plans (optional enrichment of 1.2.x)
    meth_gp = _meth_pick(meth_files, r"GroupPlan")
    gp_rows = None
    if meth_gp:
        gp_rows = read_tsv(meth_gp) or None
    meth_dmr = _meth_pick(meth_files, r"DMRPlan")
    dmr_rows = None
    if meth_dmr:
        dmr_rows = read_tsv(meth_dmr) or None

    meth_filter_html, _meth_names = section_meth_filter(meth_files, groups)

    html_parts = {}
    html_parts["s1_1"] = section_abstract(groups, individuals, total_stats)
    html_parts["s1_2_1"] = section_groups(groups, gp_rows)
    html_parts["s1_2_2"] = section_deg_comparison(args.deg_comparison, groups, dmr_rows)
    html_parts["s2"] = section_reference(args.genome, args.gtf, args.project_id,
                                         {"参考物种 (RNA)" : args.genome or ""} if args.genome else None)
    s31 = section_trimming(outdir, groups) if (outdir and groups) else na_block("samplesheet/输出目录缺失")
    if meth_filter_html:
        s31 = "<p class='cap'>RNA-seq 端过滤统计</p>" + s31 \
              + "<p class='cap'>甲基化端过滤统计 (分表)</p>" + meth_filter_html
    html_parts["s3_1"] = s31
    html_parts["s3_2"] = section_fastqc(outdir, groups) if (outdir and groups) else na_block()
    html_parts["s4_1"] = section_alignment(outdir, groups) if (outdir and groups) else na_block()
    if meth_files:
        _gm = section_meth_alignment(meth_files)
        if "暂无数据" not in _gm:
            html_parts["s4_1"] = "<p class='cap'>RNA-seq 端 (STAR)</p>" + html_parts["s4_1"] \
                                 + "<p class='cap'>甲基化端 (Bismark, 分表)</p>" + _gm
    html_parts["s4_2"] = section_csites(meth_files) if meth_files else na_block("甲基化端未提供 CsitesCoverage")
    if meth_files:
        html_parts["s5_1"] = section_meth_chr(meth_files)
        html_parts["s5_2"] = section_meth_element(meth_files)
        html_parts["s5_3"] = section_meth_number(meth_files)
    else:
        html_parts["s5_1"] = section_meth_presence(meth_files)
        html_parts["s5_2"] = na_block("未提供甲基化端结果")
        html_parts["s5_3"] = na_block("未提供甲基化端结果")
    html_parts["s6"] = section_quantification(outdir, groups) if (outdir and groups) else na_block()
    html_parts["s7_1"] = section_deseq2(args.deseq2_dir)
    html_parts["s7_2"] = section_degenerics(args.deg_table)
    html_parts["s8_1"] = section_versions(args.versions)
    html_parts["s8_2"] = section_methods(args.command_line, args.genome, args.title, args.company)
    # parent sections without own body get a short lead-in
    html_parts.setdefault("s1", "")
    html_parts.setdefault("s3", "")
    html_parts.setdefault("s4", "")
    html_parts.setdefault("s5", "")
    html_parts.setdefault("s7", "")
    html_parts.setdefault("s8", "")

    toc_items = []
    section_blocks = []
    for key, (level, title) in sections_cfg.items():
        body = html_parts.get(key, "")
        toc_items.append('<li class="l{0}"><a href="#{1}">{2}{3}</a></li>'.format(
            level, key, title, "" if body else "（无数据）"))
        if not body and not title.endswith("（无数据）"):
            pass
        section_blocks.append('<h{0} id="{1}">{2}</h{0}>{3}'.format(level, key, title, body))

    logo_html = ""
    if args.logo and os.path.isfile(args.logo):
        try:
            dst = os.path.join(args.outdir, "logo.png")
            shutil.copyfile(args.logo, dst)
            logo_html = '<img src="logo.png" alt="logo">'
        except Exception:
            pass

    html = TEMPLATE.format(
        title=esc(args.title or "RNA-seq 实验报告"),
        css=CSS,
        logo_html=logo_html,
        company=esc(args.company or ""),
        date=esc(args.date or ""),
        toc="".join(toc_items),
        sections="\n".join(section_blocks),
    )

    os.makedirs(args.outdir, exist_ok=True)
    html_path = os.path.join(args.outdir, "rnaseq_report_zh.html")
    with open(html_path, "w", encoding="utf-8") as fh:
        fh.write(html)

    pdf_path = os.path.join(args.outdir, "rnaseq_report_zh.pdf")
    ok, msg = make_pdf(html_path, pdf_path)
    print(json_summary(ok, msg, html_path, sections_cfg))
    return 0


def json_summary(pdf_ok, pdf_msg, html_path, sections_cfg):
    import json as _json
    return _json.dumps({
        "html": html_path,
        "pdf_ok": pdf_ok,
        "pdf_msg": pdf_msg,
        "sections": len(sections_cfg),
    }, ensure_ascii=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", default=".", help="pipeline output directory")
    ap.add_argument("--samplesheet", default="")
    ap.add_argument("--versions", default="")
    ap.add_argument("--deseq2-dir", default="")
    ap.add_argument("--meth-dir", default="")
    ap.add_argument("--deg-table", default="")
    ap.add_argument("--logo", default="")
    ap.add_argument("--deg-comparison", default="")
    ap.add_argument("--outdir", default="results/report_zh")
    ap.add_argument("--title", default="RNA-seq 实验报告")
    ap.add_argument("--company", default="")
    ap.add_argument("--date", default="")
    ap.add_argument("--command-line", default="")
    ap.add_argument("--genome", default="")
    ap.add_argument("--gtf", default="")
    ap.add_argument("--project-id", default="")
    ap.add_argument("--terms", default=os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "report_zh_terms.yml"))
    a = ap.parse_args()
    # normalise shell-quoted values
    for f in ("deg_table", "logo", "command_line", "title", "company", "date",
              "genome", "gtf", "project_id", "deg_comparison"):
        setattr(a, f, _unquote(getattr(a, f)))
    return render_report(a)


if __name__ == "__main__":
    sys.exit(main())
