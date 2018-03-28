rule basecall_albacore:
    input:
        "data/reads"
    output:
        "data/basecalled/workspace/pass/"
    params:
        kit=config["kit"],
        flowcell=config["flowcell"],
        save_path="data/basecalled"
    singularity:
        "containers/albacore.simg"
    threads: config["threads"]
    resources:
        mem_mb=16000
    shell:
        """
        read_fast5_basecaller.py --worker_threads {threads} \
          --save_path {params.save_path} --input {input} \
          --flowcell {params.flowcell} --kit {params.kit} --output_format fastq \
          --recursive
        """

rule concat_and_gzip_fastq:
    input:
        "data/basecalled/workspace/pass/"
    output:
        "data/basecalled/{SAMPLE}_merged.fastq.gz"
    shell:
        "cat {input}/*.fastq | gzip > {output} "
