rule demultiplex_basecall_albacore:
    input:
        "data/reads"
    output:
        "data/basecalled/workspace/pass"
    params:
        kit=config["kit"],
        flowcell=config["flowcell"],
    singularity:
        config["containers"]["albacore"]
    threads: 32
    resources:
        mem_mb=16000
    shell:
        "read_fast5_basecaller.py --barcoding --worker_threads {threads} "
        "--save_path data/basecalled --input {input} --flowcell {params.flowcell} "
        "--kit {params.kit} --output_format fastq --recursive"
