import os
from snakemake.utils import min_version

min_version("4.2.0")

configfile: "config.yaml"

# if multiplexed, add expected barcodes
MULTIPLEXED = True
if MULTIPLEXED:
    # change to appropriate barcode labels
    BARCODES = ["BC01", "BC02", "BC03", "BC04", "BC05"]
    RULES_SUBDIR = "demultiplex"
    rule all:
        input:
            expand("report_{barcode}.html", barcode=BARCODES)
else:
    # ENTER SAMPLE NAME
    SAMPLE = "test"
    RULES_SUBDIR = ""
    rule all:
        input:
            "report.html"


rules_dir = os.path.join('rules', RULES_SUBDIR)

include: os.path.join(rules_dir, 'basecall.smk')
include: os.path.join(rules_dir, 'porechop.smk')
include: os.path.join(rules_dir, 'remove_contamination.smk')
include: os.path.join(rules_dir, 'mykrobe.smk')
include: os.path.join(rules_dir, 'reports.smk')
