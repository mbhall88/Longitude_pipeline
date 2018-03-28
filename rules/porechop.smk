rule porechop:
    input:
        "data/basecalled/{sample}_merged.fastq.gz".format(sample=SAMPLE)
    output:
        "data/porechopped/{sample}_porechop.fastq.gz".format(sample=SAMPLE)
    singularity:
        "containers/nanoporeqc.simg"
    threads: config["threads"]
    log:
        "logs/porechop.log"
    shell:
        "porechop --input {input}  --output {output} --threads {threads} "
        "--check_reads 25000 --extra_end_trim 10 --discard_middle > {log}"
