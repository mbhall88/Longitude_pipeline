rule map_minimap2:
    input:
        config["tb_reference"],
        "data/porechopped/{SAMPLE}_porechop.fastq.gz"
    output:
        temp("data/mapped/{SAMPLE}.bam")
    threads:
        config["threads"]
    log:
        "logs/minimap2.log"
    singularity:
        config["containers"]["nanoporeqc"]
    resources:
        mem_mb=12000
    shell:
        "(minimap2 -t {threads} -ax map-ont {input} | "
        "samtools view -b - > {output}) 2> {log}"


rule samtools_sort:
    input:
        "data/mapped/{SAMPLE}.bam"
    output:
        "data/sorted/{SAMPLE}_sorted.bam"
    threads:
        config["threads"]
    log:
        "logs/samtools_sort.log"
    singularity:
        config["containers"]["nanoporeqc"]
    shell:
        "samtools sort -@ {threads} {input} 2> {log} > {output}"


rule samtools_index:
    input:
        "data/sorted/{SAMPLE}_sorted.bam"
    output:
        "data/sorted/{SAMPLE}_sorted.bam.bai"
    log:
        "logs/samtools_index.log"
    singularity:
        config["containers"]["nanoporeqc"]
    shell:
        "samtools index -b {input} 2> {log}"

rule bam_to_fastq:
    input:
        "data/sorted/{SAMPLE}_sorted.bam"
    output:
        "data/filtered/{SAMPLE}_filtered.fastq.gz"
    log:
        "logs/bam_to_fastq.log"
    singularity:
        config["containers"]["nanoporeqc"]
    shell:
        "samtools fastq -F 0x4 {input} > {output} 2> {log}"
