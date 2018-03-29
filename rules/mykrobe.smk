rule demultiplex_mykrobe:
    input:
        "data/filtered/{sample}_filtered.fastq.gz"
    output:
        "data/mykrobe/{sample}/{sample}_predict.json"
    params:
        species="tb",
        container=config["containers"]["mykrobe"]
    log:
        "logs/mykrobe_{sample}.log"
    shell:
        "scripts/run_mykrobe.sh {wildcards.sample} {params.species} {input} "
        "{log} {output} {params.container}"
