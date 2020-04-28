Introduction
============

CAMPAREE is an RNA expression simulator that is primed using real data to give
realistic output. For input, CAMPAREE requires fastq files and a reference
genome sequence, matching the source organism for the fastq files, with
transcript annotations. For each sample, CAMPAREE generates a simulated set of
RNA transcripts mimicking expression levels from the fastq files and accounting
for isoform-level differences and allele-specific expression. It also generates
simulated diploid genomes and their corresponding annotations with phased SNP
and indel calls from the fastq reads. Additionally CAMPAREE outputs all
underlying distributions used to simulate the expressed transcripts.

The simulated molecules are a combination of full-length RNA transcripts and
pre-mRNA transcripts with introns included. This output is saved as a text file.


About
-----

CAMPAREE is produced at the Institute for Translational Medicine and
Therapeutics, University of Pennsylvania.
