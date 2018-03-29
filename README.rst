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
**Note:** the following instructions assume you are working on a Linux operating system.

The first thing to do is download this repository onto the machine you want to run the analysis on. In the spirit of making everything reproducible and tidy I would advice to download this repository (or copy the blank version) once for each nanopore experiment you want to analyse.

Let's create our experiment directory and clone the pipeline.

.. code-block:: bash

    experiment=sample1
    git clone https://github.com/mbhall88/Longitude_pipeline.git "$experiment"

This will download the pipeline repository into a directory named, in this case, ``sample1``.

Install Singularity
---------------------
The next thing we need to do is install Singularity_. This is a program that allows for building and executing of software containers (tiny self-contained computers with pre-installed software). If you run ``singularity --help`` from the command line and get an error mesage, then you will need to install. This repository comes with a script to this, however if for some reason that doesn't work, there are more `detailed instructions here`_. **Note:** you will need to have 'root' access on the machine you are working on in order to complete the installation. If you don't have such access then contact your systems adminstrator.

.. code-block:: bash

    sudo ${experiment}/scripts/install_singularity.sh
    singularity --help

You should now get the Singularity help menu.

Install Snakemake
---------------------
Snakemake_ is a workflow management system which coordinates the running of this pipeline. In order to install it you will need to make sure you have Python3_ installed (run ``python3 -V`` to confirm this). To install, run the following.

.. code-block:: bash

    pip3 install snakemake docutils

Note: ``docutils`` is also required in order to generate the reports at the end.

Setup
========================================

.. _`Mykrobe predict`: http://www.mykrobe.com/products/predictor/
.. _Singularity: http://singularity.lbl.gov/
.. _`detailed instructions here`: http://singularity.lbl.gov/install-linux
.. _Snakemake: https://snakemake.readthedocs.io/en/stable/index.html
.. _Python3: https://www.python.org/downloads/source/
