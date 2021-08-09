.. _camparee_output:

CAMPAREE Output
===============


Using CAMPAREE with other RNA-seq simulators
--------------------------------------------

CAMPAREE is designed to simulate a realistic collection of RNA molecules. Users
can combine CAMPAREE's output with existing RNA-Seq simulators to generate FASTQ
files. This is the *in silico* equivalent of using preparing and sequencing an
RNA-Seq library from an RNA sample. Using CAMPAREE in this way can enhance
existing RNA-Seq simulators by giving them added functionality (e.g. allele-
specific expression).

Note, to perform the operations described below, CAMPAREE must be configured to
output molecules in text format.

Preparing CAMPAREE Output
^^^^^^^^^^^^^^^^^^^^^^^^^

Several existing RNA-Seq simulators, like BEERS and Polyester, allow users to
provide their own gene models and count tables to prime simulations. However,
these tools were not designed to accept CAMPAREE output directly. To solve this
problem, CAMPAREE includes the *molecule_file_to_fasta_and_count_table.py*
script to convert a CAMPAREE molecule file into a FASTA file of transcript
sequences and count table of quantification values associated with each
transcript.

This script contains several options to tailor CAMPAREE output to work with
different RNA-seq simulators. Here are the usage details for this script:

.. code-block:: none

    usage: molecule_file_to_fasta_and_count_table.py [-h] -i INPUT_MOLECULE_FILE
                                                     -o OUTPUT_DIRECTORY
                                                     [-p OUTPUT_PREFIX] [-s] [-t]

    Create transcriptome FASTA and count table from a CAMPAREE molecule file

    optional arguments:
      -h, --help            show this help message and exit
      -i INPUT_MOLECULE_FILE, --input_molecule_file INPUT_MOLECULE_FILE
                            Path to the molecule file.
      -o OUTPUT_DIRECTORY, --output_directory OUTPUT_DIRECTORY
                            Path to output directory. FASTA and count table will
                            be saved here.
      -p OUTPUT_PREFIX, --output_prefix OUTPUT_PREFIX
                            Prefix to add to all output files.
      -s, --separate_by_parent
                            Generate separate FASTA and count table files for
                            molecules from each parental genome. Assumes
                            transcript IDs in molecule file have _1/_2 suffixes to
                            identify the source parent.
      -t, --trim_polya_tails
                            Trim polyA tails from the molecules before returning
                            the FASTA. Tail length determined by soft-clipping at
                            the end of the molecule CIGAR string.

Polyester
^^^^^^^^^

Polyester can simulate RNA-Seq data given a FASTA file of transcript sequences
and a count matrix of expression values for each transcript. In order to generate
allele-specific data, users just need to provide separate transcript sequences
and quantification values for each parental copy. The instructions below
describe how to generate these files from CAMPAREE input and to generate allele-
specific RNA-Seq data with a single Polyester run, instead of running Polyester
separately for each parental allele.

1. Activate CAMPAREE environment by running the following command from the
   CAMPAREE directory::

     source ./venv_camparee/bin/activate

2. Run *molecule_file_to_fasta_and_count_table.py* script to generate a single
   FASTA file and count table containing data from both parental alleles::

     python path/to/CAMPAREE/bin/molecule_file_to_fasta_and_count_table.py \
         -i path/to/CAMPAREE/molecule_file.txt \
         -o path/to/output/dir/ \
         -p "Output_Prefix." \
         -t

3. To run Polyester, take note of the files generated in the previous step, and
   run the following R code (interactively or as a script):

     .. code-block:: R

        library(polyester)
        count_mat = read.table("path/to/CAMPAREE.count_table.txt", header = TRUE)
        simulate_experiment_countmat("path/to/CAMPAREE.transcriptome.fa",
                                     readmat=as.matrix(count_mat$Count),
                                     outdir = "path/to/Polyester/output/",
                                     seed = 100,
                                     strand_specific = TRUE,
                                     reportCoverage = TRUE)

After running, Polyester will save two FASTA files containing the simulated
RNA-Seq reads. The FASTA header lines for each read will identify their source
transcript. The "_1" and "_2" suffixes identify the parental allele for the
transcript and the "pre_mRNA" prefix identifies transcripts that retain their
intron sequences.

BEERS
^^^^^

The BEERS simulator can generate simulated allele-specific RNA-Seq data given
the following CAMPAREE output files:

- Parental genome files (FASTA)
- Parental annotation (tab-delimited text file)
- Count table (generated by  *molecule_file_to_fasta_and_count_table.py* script)

In order to generate simulated data for each allele, users need to prepare and
perform separate BEERS runs for transcripts from each parental allele.

1. Activate CAMPAREE environment by running the following command from the
   CAMPAREE directory::

     source ./venv_camparee/bin/activate

2. Run *molecule_file_to_fasta_and_count_table.py* script to generate separate
   FASTA files and count tables for each parental allele::

     python path/to/CAMPAREE/bin/molecule_file_to_fasta_and_count_table.py \
         -i path/to/CAMPAREE/molecule_file.txt \
         -o path/to/output/dir/ \
         -p "Output_Prefix." \
         -t -s

3. BEERS requires a series of configuration files that define the sequences
   and genomic locations for all simulated transcripts. It is best to create a
   set of config files for each parental genome in separate directories. The
   remaining steps assume users have created PARENTAL_GENOME_1/ and
   PARENTAL_GENOME_2/ directories::

     mkdir path/to/PARENTAL_GENOME_1/
     mkdir path/to/PARENTAL_GENOME_2/

4. CAMPAREE transcript annotation files require a few tweaks before they are
   ready for use with BEERS. Run the *prep_camparee_output_for_beers.pl*
   script on each of the parental annotation files to make these changes::

     perl path/to/prep_camparee_output_for_beers.pl \
          path/to/camparee/output/updated_annotation_1_trimmed.txt \
          > PARENTAL_GENOME_1/updated_annotation_1_edited.no_header.txt
     perl path/to/prep_camparee_output_for_beers.pl \
          path/to/camparee/output/updated_annotation_2_trimmed.txt \
          > PARENTAL_GENOME_2/updated_annotation_2_edited.no_header.txt

5. The remaining scripts use a combination of BEERS index creation scripts and
   unix commands. The commands listed here are only for one parental genome and
   need to be repeated for the second parental genome. Navigate to the directory
   of config files for the first parental genome::

     cd path/to/PARENTAL_GENOME_1/


6. Remove entries with identical coordinates and intron/exon structures (even if
   they have different gene/tx IDs). This script re-numbers the gene IDs, so
   there are no gaps left by the removed transcripts::

     perl path/to/beers/index_creation/remove_dups_in_geneinfo.pl \
          updated_annotation_1_edited.no_header.txt

7. Remove transcripts that are < 200 bp in length. This number is chosen because
   BEERS most commonly simulate Illumina PE 2x100 bp reads. The simulator
   performs best when generating reads from transcripts that are equal to or
   larger than minimum fragment length. This script re-numbers the gene IDs, so
   there are no gaps left by the removed transcripts::

     perl path/to/beers/index_creation/remove_things_too_short_in_geneinfo.pl \
          remove_dups.out

8. The fix_annotation2.pl script includes several general fixes to the gene
   annotation file that can interfere with the simulator. These include removing
   transcripts with introns < 10 bp, or transcripts with introns < 20 bp if they
   don't have canonical splice junctions. Note, this step requires a parental
   genome FASTA file prepared by CAMPAREE::

     perl path/to/beers/index_creation/fix_annotation2.pl \
          remove_things_too_short.out \
          path/to/camparee/output/custom_genome_1_edited.fa \
          > fix_annotation2.out

9. Next, since CAMPAREE also represents pre-mRNA, create separate entries for
   the pre-mRNA versions of each transcript. This allows for more granular
   control over how BEERS generates pre-mRNAs. Add these pre-mRNA sequences
   after the above filtering steps to avoid cases where a pre-mRNA sequence is
   still present, but the mature form has been filtered out (or vice versa)::

     perl path/to/add_pre_mRNA_to_annotation.pl \
          fix_annotation2.out \
          > fix_annotation2.w_pre_mRNA.out

10. Re-number BEERS gene IDs to account for anything removed during the previous
    steps::

     perl path/to/beers/index_creation/change_names_to_GENE.i_for_geneinfo.pl \
          fix_annotation2.w_pre_mRNA.out \
          > renumbered_geneids.out

11. Remove Ensembl gene and transcript IDs, leaving the generic "GENE.X"
    designators that the simulator uses. Also, map which "GENE.X" ids correspond
    to the original Ensembl IDs::

     awk 'BEGIN{OFS = "\t"; print "Simulator.ID\tCAMPAREE.Transcript\tCAMPAREE.Gene"}; {print $10, $8, $9}' \
          renumbered_geneids.out \
          > Parental_Genome_1.BEERS_trancripts_to_CAMPAREE_transcripts.txt
     cut -f 1-7,10 renumbered_geneids.out \
         > simulator_config_geneinfo_parental_genome_1_from_camparee

12. Generate a master list of exon coordinates from gene models. This script
    generates the "master_list_of_exons.txt" file used in the next step::

     perl path/to/beers/index_creation/get_master_list_of_exons_from_geneinfofile.pl \
          simulator_config_geneinfo_parental_genome_1_from_camparee

13. Generate the geneseq file from the master list of exons. The temp file
    stores an updated version of the geneinfo file if any of the gene models
    were skipped (e.g. they were not located on any chromosomes or in any region
    defined in the genome fasta file). Note: this script renames the original
    geneinfo file to match the filename specified here::

     perl path/to/beers/index_creation/make_fasta_files_for_master_list_of_genes.pl \
          path/to/camparee/output/custom_genome_1_edited.fa \
          master_list_of_exons.txt \
          simulator_config_geneinfo_parental_genome_1_from_camparee \
          temp.simulator_config_geneinfo_parental_genome_1_from_camparee \
          > simulator_config_geneseq_parental_genome_1_from_camparee
     mv temp.simulator_config_geneinfo_parental_genome_1_from_camparee \
        simulator_config_geneinfo_parental_genome_1_from_camparee

14. Generate a master list of intron coordinates from gene models. This script
    generates the "master_list_of_introns.txt" file used in the next step::

     perl path/to/beers/index_creation/get_master_list_of_introns_from_geneinfofile.pl \
          simulator_config_geneinfo_parental_genome_1_from_camparee

15. Generate the intronseq file from the master list of exons. Double-check the
    output of this file. Its error checking on the input files is not as
    stringent as some of the other scripts. For example, it doesn't check for
    the existence of the master_list_of_introns.txt file. You can give it any
    filename and if it does not exist, it generates an empty intronseq file
    without throwing any errors::

     perl path/to/beers/index_creation/make_fasta_file_for_master_list_of_introns.pl \
          path/to/camparee/output/custom_genome_1_edited.fa \
          master_list_of_introns.txt \
          > simulator_config_intronseq_parental_genome_1_from_camparee

16. Prepare gene_dist file. This file stores the probability distribution of
    expression all transcripts in the BEERS annotation. Prepare this from the
    CAMPAREE count matrix by running the following code in R:

     .. code-block:: R

        library(readr)
        library(dplyr)
        library(tidyr)
        library(tibble)

        beers_to_camparee =
            read_tsv("Parental_Genome_1.BEERS_trancripts_to_CAMPAREE_transcripts.txt",
                     col_types =
                        cols(Simulator.ID = col_character(),
                             CAMPAREE.Transcript = col_character(),
                             CAMPAREE.Gene = col_character())) %>%
            rename(Transcript_ID = CAMPAREE.Transcript,
                   Gene_ID = CAMPAREE.Gene)

        # Remove the "_1" and "_2" parental genome prefixes from the
        # transcript IDs. This way they'll match the IDs used by BEERS.
        tx_count_table =
            read_tsv("path/to/CAMPAREE.count_table.Parental_genome_1.txt",
                     col_types = cols(Transcript_ID = col_character(),
                                      Count = col_double())) %>%
            mutate(Transcript_ID = gsub("_1", "", Transcript_ID))

        # Merge the molecule counts with the BEERS transcript IDs and output the
        # gene_dist.txt file.
        beers_to_camparee_1 %>%
            left_join(tx_count_table_1,
                      by = "Transcript_ID") %>%
            select(-Transcript_ID, -Gene_ID) %>%
            mutate(Count = ifelse(is.na(Count), 0, Count)) %>%
            write_tsv(paste0("Parental_Genome_1.gene_dist.txt"),
                      col_names = FALSE, quote_escape = FALSE)

17. Generate random featurequants file from the gene distributions created
    above. Again, check the output of this script. It will run even if the
    gene_dist.txt file does not exist. In this case, it will assign 0 expression
    to all genes. So while the file may appear to be okay at first glance, it
    could be completely filled with genes that have 0 expression::

     perl path/to/beers/index_creation/make_featurequants.from_count_data.pl \
          simulator_config_geneinfo_parental_genome_1_from_camparee \
          Parental_Genome_1.gene_dist.txt \
          100 \
          > simulator_config_featurequantifications_parental_genome_1_from_camparee

18. Determine the number of reads to simulate, based upon the contents of the
    CAMPAREE count table::

     awk 'BEGIN{total=0}; (FNR > 1){total = total + $2}; END{print total}' \
         path/to/CAMPAREE.count_table.Parental_genome_1.txt

19. Having created all of the BEERS config files above, now the user needs to
    run BEERS to generate simulated reads. Again, this needs to be repeated for
    each parental genome. The BEERS command listed here is configured to run in
    an LSF cluster environment::

     perl path/to/beers/reads_simulatorP_updated3.pl \
          <NUMBER_OF_READS_TO_SIMULATE> \
          Parental_Genome_1 \
          500000 \
          -strandspecific \
          -error 0 \
          -subfreq 0 \
          -indelfreq 0 \
          -intronfreq 0 \
          -palt 0 \
          -sn \
          -configstem parental_genome_1_from_camparee \
          -customcfgdir path/to/beers/config/directory/PARENTAL_GENOME_1/ \
          -outdir path/to/output/directory/PARENTAL_GENOME_1/ \
          -fraglength 100,250,500

20. Repeat all of the above steps for other parental genomes.


RSEM
^^^^

RSEM can be used to both analyze and simulate allele-specific data, if
configured with an allele-specific index. However, it is simpler to

To run the RSEM simulator with CAMPAREE output, users can follow a similar
procedure as Polyester. That is, users generate simulated data from both
parental alleles with the FASTA file of CAMPAREE transcript sequences and the
table of transcript counts (produced by the *molecule_file_to_fasta_and_count_table.py*
script). Additionally, RSEM needs

1. Activate CAMPAREE environment by running the following command from the
   CAMPAREE directory::

     source ./venv_camparee/bin/activate

2. Run *molecule_file_to_fasta_and_count_table.py* script to generate a single
   FASTA file and count table containing data from both parental alleles::

     python path/to/CAMPAREE/bin/molecule_file_to_fasta_and_count_table.py \
         -i path/to/CAMPAREE/molecule_file.txt \
         -o path/to/output/dir/ \
         -p "Output_Prefix." \
         -t

3. Create an RSEM index from the FASTA file of parental transcript sequences.
   This step requires the user to provide an aligner. Since RSEM supports
   Bowtie2, the command here points RSEM to the Bowtie2 executable packaged in
   the CAMPAREE installation directory::

     rsem-prepare-reference \
         --bowtie2 \
         --bowtie2-path path/to/CAMPAREE/third_party_software/bowtie2-2.3.4.3-linux-x86_64 \
         path/to/CAMPAREE.transcriptome.fa \
         path/to/rsem/index/directory/camparee_output.transcriptome

4. When RSEM simulates data, it requires various models of sequencing biases
   (e.g. error profiles, fragment length distributions). While it is possible to
   create the model file defining this information from scratch, it is quite
   complex. Alternatively, it is simpler to quickly generate this model file by
   by using RSEM to quantify some of the FASTQ files used as input for CAMPAREE.
   This command uses the same RSEM transcriptome index created in the previous
   step. Note, the '--paired-end' and '--strandedness reverse' options will be
   specific to the dataset (i.e. library prep used to prepare original RNA-Seq
   libraries)::

     rsem-calculate-expression \
         --bowtie2 \
         --bowtie2-path path/to/CAMPAREE/third_party_software/bowtie2-2.3.4.3-linux-x86_64 \
         --paired-end \
         --strandedness reverse \
         path/to/input.Forward.fastq \
         path/to/input.Reverse.fastq \
         path/to/rsem/index/directory/camparee_output.transcriptome \
         path/to/rsem/output/directory/quant_from_input_fastq

   The above command creates the RSEM model file we need, in the following
   location::

     path/to/rsem/output/directory/quant_from_input_fastq.stat/quant_from_input_fastq.model

5. To read quantification data for simulations, RSEM requires the data be in a
   specific format. Specifically, the data must match the same file format RESM
   uses to output its quantification results. The following R code takes the
   transcript count table from the CAMPAREE output and reformats it to match
   the RESM \*.genes.results quantification file format. RSEM expects the
   quantification values in TPM (transcripts per million), which involves a
   normalization by transcript length. To get the length for each transcript in
   each parental allele, users will import the *updated_annotation_\?_trimmed.txt*
   files created by CAMPAREE for each parental genome. These are the same
   CAMPAREE annotation files used above when running BEERS with CAMPAREE output.

     .. code-block:: R

        library(readr)
        library(stringr)
        library(dplyr)
        library(tidyr)

        count_table_file = "path/to/CAMPAREE.count_table.txt"
        annotation_file.parent_1 = "path/to/camparee/output/updated_annotation_1_trimmed.txt"
        annotation_file.parent_2 = "path/to/camparee/output/updated_annotation_2_trimmed.txt"
        output_genes_results_file = "path/to/rsem/simulation/directory/CAMPAREE.genes.results"

        count_table =
        read_tsv(count_table_file, col_names = c("gene_id", "expected_count"),
                skip = 1)

        annotation_table.parent_1 =
        read_tsv(annotation_file.parent_1,
                col_types = cols(.default = col_character(),
                                 txStart = col_integer(),
                                 txEnd = col_integer(),
                                 exonCount = col_integer())) %>%
        mutate(gene_id = paste0(transcriptID, "_1")) %>%
        select(gene_id, txStart, txEnd, exonCount, exonStarts, exonEnds)
        annotation_table.parent_2 =
        read_tsv(annotation_file.parent_2,
                col_types = cols(.default = col_character(),
                                 txStart = col_integer(),
                                 txEnd = col_integer(),
                                 exonCount = col_integer())) %>%
        mutate(gene_id = paste0(transcriptID, "_2")) %>%
        select(gene_id, txStart, txEnd, exonCount, exonStarts, exonEnds)

        # Calculate lengths of each transcript
        annotation_table.both_parents.w_lengths =
        bind_rows(annotation_table.parent_1, annotation_table.parent_2) %>%
        # Pairs of start and stop coordinates define a one-based,
        # closed interval (i.e. [start, stop]). When subtracting
        # each start and stop, add 1 to the result to get the true
        # length of the interval in bp. For the pre_mRNA length,
        # include the introns in the length calculations.
        mutate(length.pre_mRNA = as.integer(txEnd - txStart + 1)) %>%
        rowwise() %>%
        mutate(ExonStarts.list = stringr::str_split(exonStarts, stringr::fixed(",")),
              ExonEnds.list = stringr::str_split(exonEnds, stringr::fixed(",")),
              ExonLengths = list(as.integer(ExonEnds.list) - as.integer(ExonStarts.list)),
              length.mRNA = sum(ExonLengths) + exonCount) %>%
        select(gene_id, length.pre_mRNA, length.mRNA) %>%
        arrange(gene_id) %>%
        ungroup()

        # Generate and format all information to match an rsem
        # genes.results file specifications.
        effective_len_min = 20
        camparee_derived_rsem_results =
        left_join(count_table %>%
                     mutate(gene_id.w_parent = stringr::str_remove(gene_id, "_pre_mRNA")),
                 annotation_table.both_parents.w_lengths,
                 by = c("gene_id.w_parent" = "gene_id")) %>%
        mutate(length =
                  case_when(stringr::str_detect(gene_id, "pre_mRNA") ~ length.pre_mRNA,
                            TRUE ~ length.mRNA),
              # Use 2x100bp (for paired-end reads) as fragment length
              effective_length = length - 200 + 1,
              `transcript_id(s)` = gene_id) %>%
        arrange(gene_id) %>%
        # Set genes with negative effective lengths to 0 to prevent
        # divide-by-zero errors. This is the procedure RSEM performs.
        mutate(effective_length = (effective_length > effective_len_min) * effective_length,
              expected_count = (effective_length > effective_len_min) * expected_count,
              Frac = expected_count / sum(expected_count),
              FPKM = if_else(effective_length > effective_len_min,
                             Frac * 1e9 / effective_length,
                             0),
              TPM = FPKM / sum(FPKM) * 1e6) %>%
        select(gene_id, `transcript_id(s)`, length, effective_length,
              expected_count, TPM, FPKM)

        # Save the file
        write_tsv(camparee_derived_rsem_results, output_genes_results_file)

6. Now use RSEM to generate simulated reads. This operation requires the RSEM
   index, model file, and gene.results file created in the previous steps::

     rsem-simulate-reads \
         path/to/rsem/index/directory/camparee_output.transcriptome \
         path/to/rsem/output/directory/quant_from_input_fastq.stat/quant_from_input_fastq.model \
         path/to/rsem/simulation/directory/CAMPAREE.genes.results \
         0 \
         <NUMBER_OF_READS_TO_SIMULATE> \
         path/to/rsem/simulation/directory/rsem.sim_results \
         -seed 42
