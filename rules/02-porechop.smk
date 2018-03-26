rule porechop:
    input:
        "data/basecalled/albacore_passed.fastq.gz"
    output:
        "data/porechopped/albacore_passed_porechop.fastq.gz"
    singularity:
        "containers/nanoporeqc.simg"
    threads: config["threads"]
    log:
        "logs/porechop.log"
    shell:
        "porechop --input {input}  --output {output} --threads {threads} "
        "--check_reads 25000 --extra_end_trim 10 --discard_middle > {log}"
