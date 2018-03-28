rule demultiplex_mykrobe:
    input:
        "data/filtered/{barcode}_filtered.fastq.gz"
    output:
        "data/mykrobe/{barcode}/{barcode}_predict.json"
    params:
        species="tb",
        container=config["containers"]["mykrobe"]
    log:
        "logs/mykrobe_{barcode}.log"
    shell:
        "scripts/run_mykrobe.sh {wildcards.barcode} {params.species} {input} "
        "{log} {output} {params.container}"
