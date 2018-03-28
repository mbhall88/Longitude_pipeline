rule mykrobe:
    input:
        "data/filtered/{SAMPLE}_filtered.fastq.gz"
    output:
        "data/mykrobe/{SAMPLE}_predict.json"
    params:
        species="tb",
        container=config["containers"]["mykrobe"]
    log:
        "logs/mykrobe.log"
    shell:
        "scripts/run_mykrobe.sh {wildcards.SAMPLE} {params.species} {input} "
        "{log} {output} {params.container}"
