import os
from snakemake.utils import min_version

min_version("4.2.0")

configfile: "config.yaml"

# if multiplexed, add expected barcodes
MULTIPLEXED = False
if MULTIPLEXED:
    # change to appropriate barcode labels
    BARCODES = ["BC01", "BC02", "BC03", "BC04", "BC05"]
else:
    # ENTER SAMPLE NAME
    SAMPLE = "test"

RULES_SUBDIR = ""
if MULTIPLEXED:
    rule all:
        input:
            expand("data/mykrobe/{barcode}/{barcode}_predict.json", barcode=BARCODES),
            expand("data/sorted/{barcode}_sorted.bam.bai", barcode=BARCODES)

    RULES_SUBDIR = "demultiplex"
else:
    rule all:
        input:
            "data/mykrobe/{SAMPLE}_predict.json".format(SAMPLE=SAMPLE),
            "data/sorted/{SAMPLE}_sorted.bam.bai".format(SAMPLE=SAMPLE)


rules_dir = os.path.join('rules', RULES_SUBDIR)

include: os.path.join(rules_dir, '01-basecall.smk')
include: os.path.join(rules_dir, '02-porechop.smk')
include: os.path.join(rules_dir, '03-remove_contamination.smk')
include: os.path.join(rules_dir, '04-mykrobe.smk')
