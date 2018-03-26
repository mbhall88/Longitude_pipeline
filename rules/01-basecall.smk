rule basecall_albacore:
    input:
        "data/reads"
    output:
        "data/basecalled"
    params:
        kit=config["kit"],
        flowcell=config["flowcell"],
    singularity:
        "containers/albacore.simg"
    threads: config["threads"]
    resources:
        mem_mb=16000
    shell:
        "read_fast5_basecaller.py --worker_threads {threads} "
        "--save_path {output} --input {input} --flowcell {params.flowcell} "
        "--kit {params.kit} --output_format fastq --recursive"


rule concat_and_gzip_fastq:
    input:
        "data/basecalled"
    output:
        "data/basecalled/albacore_passed.fastq.gz"
    shell:
        "cat {input}/workspace/pass/*.fastq | gzip > {output}"
