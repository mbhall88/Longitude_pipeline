import os
from snakemake.utils import min_version

min_version("4.2.0")

configfile: "config.yaml"

# if multiplexed, add expected barcodes
MULTIPLEXED = True
if MULTIPLEXED:
    # change to appropriate barcode labels
    BARCODES = ["BC01", "BC02", "BC03", "BC04", "BC05"]


RULES_SUBDIR = ""
if MULTIPLEXED:
    rule all:
        input:


    RULES_SUBDIR = "demultiplex"
else:
    rule all:
        input:
            "data/mapped/albacore_passed_porechop_sorted.bam",
            "data/mapped/albacore_passed_porechop_sorted.bam.bai"


rules_dir = os.path.join('rules', RULES_SUBDIR)

include: os.path.join(rules_dir, '01-basecall.smk')
include: os.path.join(rules_dir, '02-porechop.smk')
include: os.path.join(rules_dir, '03-remove_contamination.smk')
