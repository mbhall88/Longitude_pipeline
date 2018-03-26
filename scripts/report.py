from snakemake.utils import report

with open(snakemake.input.T1) as vcf:
    n_calls = sum(1 for l in vcf if not l.startswith("#"))

report("""
An example variant calling workflow
===================================

Reads were mapped to the Yeast
reference genome and variants were called jointly with
SAMtools/BCFtools.

This resulted in {n_calls} variants (see Table T1_).
Benchmark results for minimap2 can be found in tables T2_.
""", snakemake.output[0], **snakemake.input)
