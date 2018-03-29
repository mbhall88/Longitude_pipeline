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
    SAMPLES = barcode_parser(config["barcodes"])
else:
    SAMPLES = [config["sample_name"]]

rule all:
    input:
        expand("report_{sample}.html", sample=SAMPLES)

include: 'rules/basecall.smk'
include: 'rules/porechop.smk'
include: 'rules/remove_contamination.smk'
include: 'rules/mykrobe.smk'
include: 'rules/reports.smk'
