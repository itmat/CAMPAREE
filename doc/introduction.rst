Introduction
============

CAMPAREE is an RNA expression simulator that is primed using real data to give
realistic output. It estimates the signals as closely as possible from the real
data, which then becomes the ground truth in the simulated data. For input,
CAMPAREE requires RNA-Seq fastq files and, from the same organism, a reference
genome with transcript level annotation. For each sample, CAMPAREE generates a
simulated set of RNA transcripts mimicking expression levels from the fastq
files and accounting for isoform-level differences and allele-specific
expression. It also generates simulated diploid genomes and their corresponding
annotations with phased SNP and indel calls from the fastq reads. Additionally,
CAMPAREE outputs all underlying distributions used to simulate the expressed
transcripts. The output transcripts are full-length RNA transcripts or pre-mRNA
transcripts with introns included. If generating diploid data, since CAMPAREE
infers parental genomes, it is necessary that each sample be derived from a
single individual (rather than pooled RNA from multiple individuals).

This simulated transcriptome can then be used, for example, as input to an
RNA-Seq simulator, or to generate idealized coverage plots.

ABOUT
-----

CAMPAREE is being developed by the ITMAT Bioinformatics Laboratory in The
Institute for Translational Medicine and Therapeutics, University of Pennsylvania.
