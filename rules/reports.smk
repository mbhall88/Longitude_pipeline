rule plot_pre_filtering:
    input:
        "data/porechopped/{sample}.fastq.gz"
    output:
        "data/plots/{sample}_pre_filtering.pdf"
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/pistis_pre_filtering_{sample}.log"
    shell:
        "pistis --fastq {input} --output {output} 2> {log} "


rule plot_post_filtering:
    input:
        fastq="data/filtered/{sample}_filtered.fastq.gz",
        bam="data/sorted/{sample}_sorted.bam"
    output:
        "data/plots/{sample}_post_filtering.pdf"
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/pistis_post_filtering_{sample}.log"
    shell:
        "pistis --fastq {input.fastq} --output {output} --bam {input.bam} 2> {log} "


rule stats_pre_filtering:
    input:
        "data/porechopped/{sample}.fastq.gz"
    output:
        "data/stats/{sample}_pre_filtering.txt"
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/nanostat_pre_filtering_{sample}.log"
    threads:
        config["threads"]
    shell:
        "NanoStat --fastq {input} --name {output} --threads {threads} "
        "--readtype 1D 2> {log}"


rule stats_post_filtering:
    input:
        "data/filtered/{sample}_filtered.fastq.gz"
    output:
        "data/stats/{sample}_post_filtering.txt"
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/nanostat_post_filtering_{sample}.log"
    threads:
        config["threads"]
    shell:
        "NanoStat --fastq {input} --name {output} --threads {threads} "
        "--readtype 1D 2> {log}"


rule report:
    input:
        plot_pre="data/plots/{sample}_pre_filtering.pdf",
        plot_post="data/plots/{sample}_post_filtering.pdf",
        stats_pre="data/stats/{sample}_pre_filtering.txt",
        stats_post="data/stats/{sample}_post_filtering.txt",
        mykrobe="data/mykrobe/{sample}/{sample}_predict.json",
        porechop_log="logs/porechop.log"
    output:
        "report_{sample}.html"
    log:
        "logs/report_{sample}.log"
    params:
        sample="{sample}"
    script:
        "../scripts/report.py"
