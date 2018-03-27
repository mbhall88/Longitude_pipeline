singularity: "containers/nanoporeqc.simg"

rule demultiplex_map_minimap2:
    input:
        reads=expand("data/porechopped/{barcode}.fastq.gz", barcode=BARCODES),
        reference=config["tb_reference"]
    output:
        temp("data/mapped/{barcode}.bam")
    threads:
        config["threads"]
    log:
        "logs/minimap2_{barcode}.log"
    resources:
        mem_mb=6000
    shell:
        "(minimap2 -t {threads} -ax map-ont {input.reference} {input.reads} | "
        "samtools view -b - > {output}) 2> {log}"


rule demultiplex_samtools_sort:
    input:
        "data/mapped/{barcode}.bam"
    output:
        "data/sorted/{barcode}_sorted.bam"
    threads:
        config["threads"]
    log:
        "logs/samtools_sort_{barcode}.log"
    shell:
        "samtools sort -@ {threads} {input} 2> {log} > {output}"


rule demultiplex_samtools_index:
    input:
        "data/sorted/{barcode}_sorted.bam"
    output:
        "data/sorted/{barcode}_sorted.bam.bai"
    log:
        "logs/samtools_index_{barcode}.log"
    shell:
        "samtools index -b {input} 2> {log}"
