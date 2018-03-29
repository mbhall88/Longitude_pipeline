import os
import re
from snakemake.utils import min_version

min_version("4.2.0")

configfile: "config.yaml"


class InvalidBarcode(Exception):
    __module__ = Exception.__module__


def barcode_parser(barcodes_string):
    msg = "Barcode must be of the form BC01. That is, BC followed by 2 digits."
    regex = r'\bBC\d{2}\b'
    barcodes = barcodes_string.split()
    for barcode in barcodes:
        if not (len(barcode) == 4 and re.match(regex, barcode)):
            raise InvalidBarcode(barcode + '\n' + msg)
    return barcodes


# if multiplexed, add expected barcodes
MULTIPLEXED = config["multiplexed"]
if MULTIPLEXED:
    # change to appropriate barcode labels
    BARCODES = barcode_parser(config["barcodes"])
    RULES_SUBDIR = "demultiplex"
    rule all:
        input:
            expand("report_{barcode}.html", barcode=BARCODES)
else:
    # ENTER SAMPLE NAME
    SAMPLE = [config["sample_name"]]
    RULES_SUBDIR = ""
    rule all:
        input:
            expand("report_{sample}.html", sample=SAMPLE)


rules_dir = os.path.join('rules', RULES_SUBDIR)

include: os.path.join(rules_dir, 'basecall.smk')
include: os.path.join(rules_dir, 'porechop.smk')
include: os.path.join(rules_dir, 'remove_contamination.smk')
include: os.path.join(rules_dir, 'mykrobe.smk')
include: os.path.join(rules_dir, 'reports.smk')
