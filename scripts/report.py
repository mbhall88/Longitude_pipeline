import json
import os
from snakemake.utils import report


def mykrobe_overview(filepath: str) -> str:
    """Extracts the susceptiblity information from mykrobe predict json.

    :param filepath: path to mykrobe predict output json file.

    :returns A json formatted string showing drug name and susceptibility. If
    susceptibility is resistant 'R', then information on the variant is given.

    """
    sample_id = os.path.basename(filepath).split('.')[0].split('_')[0]
    with open(filepath, 'r') as mykrobe_json:
        data = json.load(mykrobe_json)
    return json.dumps(data[sample_id]['susceptibility'], indent=4)


def get_num_reads(stats_file: str) -> int:
    """Extracts the number of reads from a nanostats text file."""
    with open(stats_file, 'r') as stats:
        for line in stats:
            if 'Number of reads:' in line:
                num_reads = line.split()[-1]
    return int(num_reads)


mykrobe_report = mykrobe_overview(snakemake.input.mykrobe)
num_reads_pre_filter = get_num_reads(snakemake.input.stats_pre)
num_reads_post_filter = get_num_reads(snakemake.input.stats_post)
percent_reads_unmapped = num_reads_pre_filter / num_reads_post_filter * 100

with open(snakemake.input.porechop_log, 'r') as log_file:
    porechop_log = log_file.read()


report("""
Report for {sample}
===================================

1. ``porechop`` was run to trim adapter sequences from reads and discard reads
with adapters found in the middle. More detailing information can be found in
porechop_log_.


""", snakemake.output, sample=snakemake.params.sample, **snakemake.input)
