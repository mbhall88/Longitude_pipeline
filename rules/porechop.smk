rule porechop:
    input:
        "data/basecalled/{sample}.fastq.gz"
    output:
        "data/porechopped/{sample}.fastq.gz"
    singularity:
        "containers/nanoporeqc.simg"
    threads: config["threads"]
    params:
        output_type=("--barcode_dir data/porechopped" if MULTIPLEXED
                     else "--output {output}"),
        unassigned=("--discard_unassigned" if MULTIPLEXED else "")
    log:
        "logs/porechop.log"
    shell:
        "porechop --input {input}  {params.output_type} --threads {threads} "
        "--check_reads 25000 --extra_end_trim 10 --discard_middle "
        "{params.unassigned} --format fastq.gz > {log}"
