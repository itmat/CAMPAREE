.. _quick-start-guide:

Quick Start Guide
=================

This guide will walk you through basic installation and usage of CAMPAREE, by
running a simulation on a simplified dataset consisting of a mouse genome
truncated to about 6 million bases and two samples of reads that align there.
You can find more detailed installation instructions :ref:`here <installing-general>`.

Installation
------------

Make sure you have the following installed on your system:

- python version 3.11
- Java 1.8

1. Create a Python virtual environment to install required Python libraries to:

    python3 -m venv ./venv_camparee

2. Activate the environment:

    source ./venv_camparee/bin/activate

3. Install CAMPAREE:

    pip install git+https://github.com/itmat/CAMPAREE

4. [Optional] Validate installation using the 'Baby Genome' as described in the :ref:`quick start guide <quick-start-baby-genome>`.

5. Prepare/install resource files for the organism of choice (see full instructions in :ref:`Resource Files <resource-pre-built>` section).

.. _quick-start-baby-genome:

Baby Genome
-----------

The "baby genome" is a truncated version of mm10 consisting of segments of length at most 1 million bases chosen from chromosomes 1, 2, 3, X, Y, and MT.

CAMPAREE comes with scripts to work with this.
First, we create the files and then we have to make a STAR index for it, using built-in scripts available when the CAMPAREE environment is activate.

    create_camparee_baby_config
    create_star_index_for_baby_genome

Perform Test Run
^^^^^^^^^^^^^^^^

We are now ready to run CAMPAREE on a two small sample fastq files aligning to the baby genome.

The default config file for the baby genome has CAMPAREE run all operations serially on a single machine.
To perform the test run with these defaults, first download the config file and then run CAMPAREE with it::

    camparee --config baby.config.yaml --run_id 1

The argument `--run_id 1` indicates that the run number is 1.
If you run this again, you must either remove the output directory `test_data/results/run_1/` or specify a new run number.

It is also possible to test deployment to a cluster.

For LSF clusters run::

    camparee --config baby.config.yaml --run_id 1 --scheduler_mode lsf

For SGE clusters run::

    camparee --config baby.config.yaml --run_id 1 --scheduler_mode sge

Check Results
^^^^^^^^^^^^^

When the run completes, all output will be saved to
``CAMPAREE/test_data/results/run_1/``. The final outputs will be in the text
files ``test_data/results/run_1/CAMPAREE/data/sample1/molecule_file.txt`` and
``test_data/results/run_1/CAMPAREE/data/sample2/molecule_file.txt``. Each line
(after the header line) corresponds to the sequence of a single simulated
molecule in a tab-delimited format. The default config file outputs 10000
molecules.

Validate the contents of these files by running the following md5sum command in
the CAMPAREE install directory::

    md5sum test_data/results/run_1/CAMPAREE/data/sample*/molecule_file.txt

The output of this command should match the following md5sum values::

    60ab3c4fd245ad3003a8113d8edbdb2b  test_data/results/run_1/CAMPAREE/data/sample1/molecule_file.txt
    6a5a1a2d9632f5fb3f3ea7740cee7514  test_data/results/run_1/CAMPAREE/data/sample2/molecule_file.txt
