rule plot_pre_filtering:
    input:
        expand("data/porechopped/{barcode}.fastq.gz", barcode=BARCODES)
    output:
        "data/plots/{barcode}_pre_filtering.pdf"
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/pistis_pre_filtering_{barcode}.log"
    shell:
        "pistis --fastq {input} --output {output} 2> {log} "


rule plot_post_filtering:
    input:
        fastq="data/filtered/{barcode}_filtered.fastq.gz",
        bam="data/sorted/{barcode}_sorted.bam"
    output:
        "data/plots/{barcode}_post_filtering.pdf"
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/pistis_post_filtering_{barcode}.log"
    shell:
        "pistis --fastq {input.fastq} --output {output} --bam {input.bam} 2> {log} "


rule stats_pre_filtering:
    input:
        expand("data/porechopped/{barcode}.fastq.gz", barcode=BARCODES)
    output:
        "data/stats/{barcode}_pre_filtering.txt"
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/nanostat_pre_filtering_{barcode}.log"
    threads:
        config["threads"]
    shell:
        "NanoStat --fastq {input} --name {output} --threads {threads} "
        "--readtype 1D 2> {log}"


rule stats_post_filtering:
    input:
        "data/filtered/{barcode}_filtered.fastq.gz"
    output:
        "data/stats/{barcode}_post_filtering.txt"
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/nanostat_post_filtering_{barcode}.log"
    threads:
        config["threads"]
    shell:
        "NanoStat --fastq {input} --name {output} --threads {threads} "
        "--readtype 1D 2> {log}"


rule report:
    input:
        plot_pre="data/plots/{barcode}_pre_filtering.pdf",
        plot_post="data/plots/{barcode}_post_filtering.pdf",
        stats_pre="data/stats/{barcode}_pre_filtering.txt",
        stats_post="data/stats/{barcode}_post_filtering.txt",
        mykrobe="data/mykrobe/{barcode}_predict.json",
        porechop_log="logs/porechop.log"
    output:
        "report_{barcode}.html"
    params:
        sample="{barcode}"
    script:
        "../scripts/report.py"
