rule demultiplex_porechop:
    input:
        "data/basecalled/workspace/pass"
    output:
        expand("data/porechopped/{barcode}.fastq.gz", barcode=BARCODES)
    singularity:
        "containers/nanoporeqc.simg"
    threads: config["threads"]
    log:
        "logs/porechop.log"
    shell:
        "porechop --input {input}  --barcode_dir data/porechopped "
        "--threads {threads} --discard_unassigned --check_reads 25000 "
        "--extra_end_trim 10 --discard_middle --format fastq.gz > {log}"
