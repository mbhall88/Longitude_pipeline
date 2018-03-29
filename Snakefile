import os
import re
from typing import List
from snakemake.utils import min_version


# need this version as minimum as it is when singularity support was added
min_version("4.2.0")


configfile: "config.yaml"


class InvalidBarcode(Exception):
    __module__ = Exception.__module__


def barcode_parser(barcodes_string: str) -> List[str]:
    """Parses the barcodes string and ensures they follow correct format"""
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


# only run basecalling when requested
if config["basecall"]:
    include: 'rules/basecall.smk'


# the snakemake files that run the different parts of the pipeline
include: 'rules/porechop.smk'
include: 'rules/remove_contamination.smk'
include: 'rules/mykrobe.smk'
include: 'rules/reports.smk'
