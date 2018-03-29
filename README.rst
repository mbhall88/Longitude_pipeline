========================================
Analysis pipeline for *M. tuberculosis*
========================================
This repository is designed to analyse Oxford Nanopore Technologies sequence data.
Specifically, for drug resistance prediction with `Mykrobe predict`_.

To achieve this, the pipeline does the following:

1. Basecalling of raw nanopore sequencing files (if required).
2. Adapter trimming of the basecalled reads (and demultiplexing if required).
3. Alignment to the *M. tuberculosis* reference genome and removal of unmapped reads.
4. Drug resistance prediction with Mykrobe predict.
5. Final report with plots and statistics, along with the Mykrobe results.

Installation
========================================
The first thing to do is download this repository onto the machine you want to run the analysis on. In the spirit of making everything reproducible and tidy I would advice to download this repository (or copy the blank version) once for each nanopore experiment you want to analyse.

Let's create our experiment directory and clone the pipeline.

.. code-block:: bash

    experiment=sample1
    git clone https://github.com/mbhall88/Longitude_pipeline.git "$experiment"

This will download the pipeline repository into a directory named, in this case, ``sample1``.



.. _`Mykrobe predict`: http://www.mykrobe.com/products/predictor/
