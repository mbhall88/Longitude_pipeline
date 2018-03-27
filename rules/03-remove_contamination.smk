singularity: "containers/nanoporeqc.simg"

rule map_minimap2:
    input:
        reads="data/porechopped/albacore_passed_porechop.fastq.gz",
        reference=config["tb_reference"]
    output:
        temp("data/mapped/albacore_passed_porechop.bam")
    threads:
        config["threads"]
    log:
        "logs/minimap2.log"
    resources:
        mem_mb=16000
    shell:
        "minimap2 -t {threads} -ax map-ont {input.reference} {input.reads} "
        "2> {log} | samtools view -b - > {output}"


rule samtools_sort:
    input:
        "data/mapped/albacore_passed_porechop.bam"
    output:
        "data/mapped/albacore_passed_porechop_sorted.bam"
    threads:
        config["threads"]
    log:
        "logs/samtools_sort.log"
    shell:
        "samtools sort -@ {threads} {input} 2> {log} > {output}"


rule samtools_index:
    input:
        "data/mapped/albacore_passed_porechop_sorted.bam"
    output:
        "data/mapped/albacore_passed_porechop_sorted.bam.bai"
    log:
        "logs/samtools_index.log"
    shell:
        "samtools index -b {input} 2> {log}"
