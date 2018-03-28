rule plot_pre_filtering:
    input:
        "data/porechopped/{sample}_porechop.fastq.gz".format(sample=SAMPLE)
    output:
        "data/plots/{sample}_pre_filtering.pdf".format(sample=SAMPLE)
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/pistis_pre_filtering.log"
    shell:
        "pistis --fastq {input} --output {output} 2> {log} "


rule plot_post_filtering:
    input:
        fastq="data/filtered/{sample}_filtered.fastq.gz".format(sample=SAMPLE),
        bam="data/sorted/{sample}_sorted.bam".format(sample=SAMPLE)
    output:
        "data/plots/{sample}_post_filtering.pdf".format(sample=SAMPLE)
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/pistis_post_filtering.log"
    shell:
        "pistis --fastq {input.fastq} --output {output} --bam {input.bam} 2> {log} "


rule stats_pre_filtering:
    input:
        "data/porechopped/{sample}_porechop.fastq.gz".format(sample=SAMPLE)
    output:
        "data/stats/{sample}_pre_filtering.txt".format(sample=SAMPLE)
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/nanostat_pre_filtering.log"
    threads:
        config["threads"]
    shell:
        "NanoStat --fastq {input} --name {output} --threads {threads} "
        "--readtype 1D 2> {log}"


rule stats_post_filtering:
    input:
        "data/filtered/{sample}_filtered.fastq.gz".format(sample=SAMPLE)
    output:
        "data/stats/{sample}_post_filtering.txt".format(sample=SAMPLE)
    singularity:
        config["containers"]["nanoporeqc"]
    log:
        "logs/nanostat_post_filtering.log"
    threads:
        config["threads"]
    shell:
        "NanoStat --fastq {input} --name {output} --threads {threads} "
        "--readtype 1D 2> {log}"


rule report:
    input:
        plot_pre="data/plots/{sample}_pre_filtering.pdf".format(sample=SAMPLE),
        plot_post="data/plots/{sample}_post_filtering.pdf".format(sample=SAMPLE),
        stats_pre="data/stats/{sample}_pre_filtering.txt".format(sample=SAMPLE),
        stats_post="data/stats/{sample}_post_filtering.txt".format(sample=SAMPLE),
        mykrobe="data/mykrobe/{sample}_predict.json".format(sample=SAMPLE),
        porechop_log="logs/porechop.log"
    output:
        "report.html"
    params:
        sample=SAMPLE
    script:
        "../scripts/report.py"
