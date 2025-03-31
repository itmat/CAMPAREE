# CAMPAREE

## Introduction

CAMPAREE is a RNA expression simulator that is primed using real data to give realistic output.
CAMPAREE needs as input a reference genome with transcript annotations as well as fastq files of samples of the species to base the output on.
For each sample, CAMPAREE outputs a simulated set of RNA transcripts mimicking expression levels with in the fastq files and accounting for isoform-level expression and allele-specific expression.
It also outputs simulated diploid genomes and their corresponding annotations with phased SNP and indel calls in the transcriptome from fastq reads.
Additionally the simulation outputs the underlying distributions used for expressing the transcripts.

## Documentation

Our documentation is on [ReadTheDocs](https://camparee.readthedocs.io).

## Quick Start Guide

This guide will walk you through basic installation and usage of CAMPAREE running a simulation on a simplified dataset consisting of a mouse genome truncated to about 6 million bases and two samples of reads that align there.

### Installation

Make sure you have the following installed on your system:

- python version 3.11
- Java 1.8

Create a Python virtual environment to install required Python libraries to:

    python3 -m venv ./venv_camparee

And activate the environment:

    source ./venv_camparee/bin/activate

Install required libraries:

    pip install git+https://github.com/itmat/CAMPAREE

### Baby Genome

The "baby genome" is a truncated version of mm10 consisting of segments of length at most 1 million bases chosen from chromosomes 1, 2, 3, X, Y, and MT.

CAMPAREE comes with scripts to work with this, but first we have to make a STAR index for it, using a built-in script available when the CAMPAREE environment is activate.

    create_star_index_for_baby_genome

### Perform Test Run

We are now ready to run CAMPAREE on a two small sample fastq files aligning to the baby genome.
If you have not already done so for installation, activate the python environment:

    source ./venv_camparee/bin/activate

The default config file for the baby genome has CAMPAREE run all operations serially on a single machine.
To perform the test run with these defaults, first download the config file and then run CAMPAREE with it:

    create_camparee_baby_config
    camparee --config baby.config.yaml --run_id 1

The argument `--run_id 1` indicates that the run number is 1.
If you run this again, you must either remove the output directory `test_data/results/run_1/` or specify a new run number.

It is also possible to test deployment to a cluster.
For LSF clusters run:

    camparee --config baby.config.yaml --run_id 1 --scheduler_mode lsf

For SGE clusters run:

    camparee --config baby.config.yaml --run_id 1 --scheduler_mode sge

### Check Results

When the run completes, output will be created in `test_data/results/run_1/`.
The final outputs will be in the text files `test_data/results/run_1/CAMPAREE/data/sample1/molecule_file` and  `test_data/results/run_1/CAMPAREE/data/sample2/molecule_file`.
Each line (after the header line) corresponds to a sequence of a single molecule in a tab-separated format.
The default config file outputs 10000 molecules.

## Packaged Software

CAMPAREE includes (in the 'third_party_software' directory) binary executables from the following pieces of software. Please cite the listed papers when using CAMPAREE.

- STAR
    * Version: 2.5.2a
    * Source code: https://github.com/alexdobin/STAR/tree/2.5.2a
    * Citation: Dobin A, Davis CA, Schlesinger F, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*. 2013 Jan 1;29(1):15-21. doi: 10.1093/bioinformatics/bts635
- BEAGLE
    * Version: 5.0 (28Sep18.793)
    * Source code: https://faculty.washington.edu/browning/beagle/beagle.html
    * Citation:
        Browning SR, Browning BL. Rapid and accurate haplotype phasing and missing data inference for whole genome association studies by use of localized haplotype clustering. *Am J Hum Genet*. 2007 Nov;81(5):1084-97. doi:10.1086/521987
- Kallisto
    * Version: 0.45.0
    * Source code: https://github.com/pachterlab/kallisto/tree/v0.45.0
    * Citation: Bray NL, Pimentel H, Melsted P, Pachter L. Near-optimal probabilistic RNA-seq quantification. *Nat Biotechnol*. 2016 May;34(5):525-7. doi: 10.1038/nbt.3519
- Bowtie2
    * Version: 2.3.4.3
    * Source code: https://github.com/BenLangmead/bowtie2/tree/v2.3.4.3
    * Citation: Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. *Nat Methods*. 2012 Mar 4;9(4):357-9. doi: 10.1038/nmeth.1923

## Funding

Work on CAMPAREE is supported by R21-LM012763-01A1: “The Next Generation of RNA-Seq Simulators for Benchmarking Analyses” (PI: Gregory R. Grant).
