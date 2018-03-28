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

reference = os.path.splitext(os.path.basename(snakemake.config["tb_reference"]))[0]
mykrobe_report = mykrobe_overview(snakemake.input.mykrobe)
num_reads_pre_filter = get_num_reads(snakemake.input.stats_pre)
num_reads_post_filter = get_num_reads(snakemake.input.stats_post)
percent_reads_mapped = round(num_reads_post_filter / num_reads_pre_filter * 100, 2)
sample=snakemake.params.sample

with open(snakemake.input.porechop_log, 'r') as log_file:
    porechop_log = log_file.read()


report("""
===================================
Report for {sample}
===================================

Quality Control
===================================
1. Porechop_ was run to trim adapter sequences from reads and discard reads with adapters found in the middle. More detailing information can be found in `porechop_log`_. For quality control plots of the reads after this step, see `plot_pre`_.
2. Reads were aligned to the TB reference {reference} using Minimap2_.
3. All reads which did not map to {reference} were removed. Prior to filtering there were {num_reads_pre_filter} reads. After filtering there remains {num_reads_post_filter}. This means {percent_reads_mapped}% of reads mapped to {reference}. For more stats on the pre-filtered reads see `stats_pre`_ and for post-filtered reads see `stats_post`_. For quality control plots of the reads after this step (and read percent identity to {reference}) see `plot_post`_. Stats were produced with NanoStat_ and plots with Pistis_.

Mykrobe Analysis
===================================
A summary of the susceptiblity information from `Mykrobe predict`_ is shown here. For the full report, see mykrobe_. 'S' means susceptible and 'R' means resistant. If resistance is identified for a drug then the predicted responsible variant is given, along with supporting information.

{mykrobe_report}


.. _Porechop: https://github.com/rrwick/Porechop
.. _Minimap2: https://github.com/lh3/minimap2
.. _NanoStat: https://github.com/wdecoster/nanostat
.. _Pistis: https://github.com/mbhall88/pistis
.. _`Mykrobe predict`: http://www.mykrobe.com/products/predictor/
""", snakemake.output[0], metadata="Author: Michael Hall (michael.hall@ebi.ac.uk)",
**snakemake.input)
