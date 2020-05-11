Resource Files
==============


Overview
--------

In order to simulate a realistic transcriptome, CAMPAREE requires information
about an organism's genome sequence, transcript structures, and karyotype.
CAMPAREE reads this information from several user-provided "resource files." By
default, these files are stored in the ``resources/`` subdirectory of the
CAMPAREE installation directory, though this can be changed in the CAMPAREE
config file. A collection of resource files is derived from a specific build of
an organism's genome and corresponding transcript annotation. We have prepared
collections of pre-built resource files for several commonly-sequenced
organisms, available for download and use with CAMPAREE. Additionally, CAMPAREE
comes bundled with utilities to generate custom resource files for users who
want to simulated data from organisms/genome builds not available as pre-built
resource files.


File Descriptions
-----------------

For organizational purposes, all CAMPAREE resource files should be associated
with a specific species/model identifier, that ideally includes the organism
name, as well as the genome and annotations builds and versions. For example,
the pre-built resource files for *MusMusculus_GRCm38_Ensemblv99* were generated
from the GRCm38 build of the mouse genome sequence and corresponding transcript
models from version 99 of Ensembl. Here are specifics on the formatting for each
of the different resource files.

Reference Genome Sequence File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The reference genome sequence is stored as plain-text in FASTA format. Each
entry in the file corresponds to a single chromosome or contig. The sequence for
each chromosome/contig is stored on a single line. The header line (starting
with a ">" character) for each entry contains the chromosome/contig identifier
and must match the chromosome/contig identifiers used in both the annotation and
chromosome ploidy files. CAMPAREE can process this file if it is uncompressed or
if it is gzipped. Gzipped files must end with the ".gz" extension to be
processed correctly.

Transcript Annotation File
^^^^^^^^^^^^^^^^^^^^^^^^^^

Gene/transcript models are stored as plain-text in a CAMPAREE-specific,
tab-delimited format, derived from the UCSC BED format. The annotation file
contains eleven columns:

1. chromosome/contig identifier
    match identifiers used in sequence and ploidy files
2. strand
    Genomic strand from which the transcript is transcribed. Either "+" or "-".
3. txStart
    Start coordinate of transcript. 1-based integer identifying the position of
    the transcript's first nucleotide in the chromosome/contig sequence stored
    in the genome sequence file.
4. txEnd
    End coordinate of transcript. 1-based integer identifying the position of
    the transcript's last nucleotide in the chromosome/contig sequence.
5. exonCount
    Number of exons in the transcript.
6. exonStarts
    Comma-delimited list of start coordinates for each exon in the transcript.
    Each coordinate is a 1-based integer identifying the position an exon's
    first nucleotide in the chromosome/contig sequence.
7. exonEnds
    Comma-delimited list of end coordinates for each exon in the transcript.
    Each coordinate is a 1-based integer identifying the position an exon's
    last nucleotide in the chromosome/contig sequence.
8. transcript_id
    Identifier for the current transcript. This must be unique across all entries
    in the annotation file. An example would be an Ensembl transcript id:
    ENST00000389707.
9. gene_id
    Identifier for the gene corresponding to the current transcript. Multiple
    transcript's can correspond to the same gene (e.g. alternative splicing),
    so gene identifiers do not need to be unique.  An example would be an Ensembl
    gene id: ENSG00000133794.
10. genesymbol
     Common gene symbol for the current transcript, or 'None'. This is used to
     facilitate downstream interpretation and use of the data, as the gene
     symbols are more commonly known than transcript / gene ids. The values in
     this column are not used as part of any processing performed by CAMPAREE.
11. biotype
     Biotype / functional classification for current gene, or 'None'. This is
     not currently used by CAMPAREE.

Chromosome Ploidy File
^^^^^^^^^^^^^^^^^^^^^^

The chromosome ploidy file is a tab-delimited text file defining the number of
copies of each chromosomes present in each gender of the simulated organism's
species. Since CAMPAREE simulates transcripts arising from a *diploid* genome,
it needs the organization and number of chromosomes in each gender of the
simulated organisms. While most chromosomes will have two copies, regardless of
gender, the sex and mitochondrial chromosomes often deviate from this pattern in
a species-dependent fashion. The chromosome/contig identifiers in this file must
match those used in the genome sequence and annotation files. However, any
chromosome/contig identifiers omitted from the ploidy file will not be used by
CAMPAREE (i.e. no transcripts will be simulated from omitted
chromosome/contigs). This allows users broad control over which
chromosomes/contigs CAMPAREE uses. For example, all of the pre-built resource
files use ploidy files that limit CAMPAREE to the standard chromosomes
(numbered/lettered, sex, and mitochondrial). For examples, please refer to the
ploidy file for the baby genome, packaged in the ``resources/baby_genome.mm10/``
subdirectory of the CAMPAREE install directory, or to one of the ploidy files
packaged with the pre-built resource files.

STAR Genome Index
^^^^^^^^^^^^^^^^^


Pre-Built Resource Files
------------------------

We have prepared resource files for the genomes of several commonly studied
organisms. The STAR Genome Indexes occupy substantially more disk space than the
other files, so we packaged them separately. The resource files are also packaged
with a README file, which contains all genome/annotation version information a
user would need to provide when reporting CAMPAREE results in a publication.
This includes md5sums for each resource file, the command used to build the
resource files, the command use to build the STAR Genome Index, and the sources/
URLs and versions for the reference genome and gene models used to build the
resource files. Effectively, the README file contains all information a user

To use the pre-build resource files, download and unpackage them in the
``resources/`` subdirectory of the CAMPAREE installation directory, and then
upodate the *resources* section of the config file to point to each of the
resource files. Here are example commands for download and installing resource
files for the mouse genome with Ensembl gene models:

.. code-block:: none

    # Navigate to CAMPAREE resources directory
    cd /path/to/CAMPAREE/resources

    # Download files
    wget https://itmat.data-simulators.s3.amazonaws.com/BEERS2/CAMPAREE_RESOURCE_FILES/MusMusculus_GRCm38_Ensemblv99__Resource_files.tar.gz ./
    wget https://itmat.data-simulators.s3.amazonaws.com/BEERS2/CAMPAREE_RESOURCE_FILES/MusMusculus_GRCm38_Ensemblv99__STAR_index.tar.gz ./

    # Unpack resource files and STAR index
    tar -xvzf MusMusculus_GRCm38_Ensemblv99__Resource_files.tar.gz
    tar -xvzf MusMusculus_GRCm38_Ensemblv99__STAR_index.tar.gz


Download Links for Pre-Built Resource Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

MusMusculus_GRCm38_Ensemblv99 (Built 2020-03-30)

- `Resource files <https://itmat.data-simulators.s3.amazonaws.com/BEERS2/CAMPAREE_RESOURCE_FILES/MusMusculus_GRCm38_Ensemblv99__Resource_files.tar.gz>`_
- `STAR index <https://itmat.data-simulators.s3.amazonaws.com/BEERS2/CAMPAREE_RESOURCE_FILES/MusMusculus_GRCm38_Ensemblv99__STAR_index.tar.gz>`_

HomoSapiens_GRCh38_Ensemblv99 (Built 2020-03-30)

- `Resource files <https://itmat.data-simulators.s3.amazonaws.com/BEERS2/CAMPAREE_RESOURCE_FILES/HomoSapiens_GRCh38_Ensemblv99__Resource_files.tar.gz>`_
- `STAR index <https://itmat.data-simulators.s3.amazonaws.com/BEERS2/CAMPAREE_RESOURCE_FILES/HomoSapiens_GRCh38_Ensemblv99__STAR_index.tar.gz>`_


Generating Custom Resource Files
--------------------------------
