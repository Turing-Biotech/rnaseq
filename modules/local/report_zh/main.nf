process REPORT_ZH {
    tag "report_zh"
    label "process_single"

    conda "${moduleDir}/environment.yml"
    container "community.wave.seqera.io/library/python:3.14.5--dc8358b3c5eeb927"

    input:
    path samplesheet
    path versions_yml
    path deseq2_files, optional: true
    path meth_files, optional: true
    path deg_table, optional: true
    path logo, optional: true
    val outdir
    val meta

    output:
    path "results/report_zh/rnaseq_report_zh.html", emit: html
    path "results/report_zh/rnaseq_report_zh.pdf", emit: pdf, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def cl = (meta.command_line ?: '').toString().replace("'", "'\\''")
    def deg = deg_table ? deg_table.toString() : ''
    def logopath = logo ? logo.toString() : ''
    '''
    set -e
    mkdir -p results/report_zh
    if [ "$logopath" != "" ] && [ -f "$logopath" ]; then
        cp -f "$logopath" results/report_zh/logo.png
    fi

    # Dependencies: pure-python wheels (jinja2, pyyaml) and best-effort PDF stack
    python -m pip install --quiet --disable-pip-version-check jinja2 pyyaml
    if ! command -v wkhtmltopdf > /dev/null 2>&1; then
        ( apt-get update -qq && apt-get install -y -qq --no-install-recommends wkhtmltopdf fonts-noto-cjk ) > /dev/null 2>&1 || true
    fi

    python ${moduleDir}/bin/render_report_zh.py \\
        --results "$outdir" \\
        --samplesheet "$samplesheet" \\
        --versions "$versions_yml" \\
        --deseq2-dir "." \\
        --meth-dir "meth" \\
        --deg-table "'$deg'" \\
        --logo "'$logopath'" \\
        --outdir "results/report_zh" \\
        --title "'${meta.title ?: 'RNA-seq 实验报告'}'" \\
        --company "'${meta.company ?: ''}'" \\
        --date "'${meta.date ?: ''}'" \\
        --command-line "'$cl'" \\
        --genome "'${meta.genome ?: ''}'" \\
        --gtf "'${meta.gtf ?: ''}'" \\
        --project-id "'${meta.project_id ?: ''}'"

    printf 'report_zh:\n    version: 1.0.0\n' > "$work/versions.yml"
    '''

    stub:
    '''
    mkdir -p results/report_zh
    touch results/report_zh/rnaseq_report_zh.html
    printf 'report_zh:\n    version: 1.0.0\n' > "$work/versions.yml"
    '''
}
