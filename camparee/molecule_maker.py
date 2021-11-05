import os
import sys
import collections
import argparse
import numpy
import pickle

from camparee.abstract_camparee_step import AbstractCampareeStep
from camparee.camparee_constants import CAMPAREE_CONSTANTS

from beers_utils.molecule_packet import MoleculePacket
from beers_utils.molecule import Molecule
from beers_utils.sample import Sample
from beers_utils.general_utils import GeneralUtils

class MoleculeMakerStep(AbstractCampareeStep):
    """
    MoleculeMaker generates molecules based off of gene, intron, and allelic
    quantification files as well as customized genomic sequence and annotation
    """

    OUTPUT_OPTIONS_W_EXTENSIONS=CAMPAREE_CONSTANTS.MOLECULE_MAKER_OUTPUT_OPTIONS_W_EXTENSIONS
    OUTPUT_FILENAME_PATTERN=CAMPAREE_CONSTANTS.MOLECULE_MAKER_OUTPUT_FILENAME_PATTERN
    DEFAULT_MOLECULES_PER_PACKET=CAMPAREE_CONSTANTS.MOLECULE_MAKER_DEFAULT_NUM_MOLECULES_PER_PACKET

    # Default filename patterns for output from previous CAMPAREE steps.
    # _GENE_QUANT_FILENAME=CAMPAREE_CONSTANTS.TXQUANT_OUTPUT_GENE_FILENAME
    # _INTRON_QUANT_FILENAME=CAMPAREE_CONSTANTS.INTRON_OUTPUT_FILENAME
    # _TX_QUANT_PSI_FILENAME=CAMPAREE_CONSTANTS.TXQUANT_OUTPUT_PSI_FILENAME
    #_ALLELIC_IMBALANCE_FILENAME=CAMPAREE_CONSTANTS.ALLELIC_IMBALANCE_OUTPUT_FILENAME
    _PARENTAL_TX_FASTA_FILENAME_PATTERN=CAMPAREE_CONSTANTS.TRANSCRIPTOME_FASTA_OUTPUT_FILENAME_PATTERN
    _PARENTAL_ANNOT_FILENAME_PATTERN=CAMPAREE_CONSTANTS.UPDATEANNOT_OUTPUT_FILENAME_PATTERN
    _PARENTAL_GENOME_FASTA_FILENAME_PATTERN=CAMPAREE_CONSTANTS.GENOMEBUILDER_SEQUENCE_FILENAME_PATTERN
    _PARENTAL_GENOME_INDEL_FILENAME_PATTERN=CAMPAREE_CONSTANTS.GENOMEBUILDER_INDEL_FILENAME_PATTERN

    def __init__(self, log_directory_path, data_directory_path, parameters=None):
        """Constructor for MoleculeMakerStep object.

        Parameters
        ----------
        data_directory_path: string
            Full path to data directory
        log_directory_path : string
            Full path to log directory.
        parameters : dict
            Dictionary of other parameters specified by the config file. This
            parameter is not used by this class and is retained for uniformity
            with all other CAMPAREE steps.

        """
        self.data_directory_path = data_directory_path
        self.log_directory_path = log_directory_path

    # Nearly all of the validation for this step is already performed in the
    # expression_pipeline, since the output options in the config file are
    # specified outside of the standard 'steps' framework.
    def validate(self):
        return True

    def load_annotation(self, file_path):
        transcripts = dict()
        with open(file_path) as annotation_file:
            for line in annotation_file:
                if line.startswith("#"):
                    continue # Comment/header line

                chrom, strand, tx_start, tx_end, exon_count, exon_starts, exon_ends, transcript_id, gene_id, gene_sybmol, *other \
                        = line.split("\t")

                transcripts[transcript_id] = (chrom, strand, int(tx_start), int(tx_end),
                                                [int(start) for start in exon_starts.split(",")],
                                                [int(end) for end in exon_ends.split(",")])
        return transcripts


    def load_transcriptome(self, file_path):
        """
        Read in a fasta file and load is a dictionary id -> sequence
        assumed one-line for the whole contig
        """
        transcripts = dict()
        with open(file_path) as transcriptome_file:
            while True:
                line = transcriptome_file.readline()
                if not line:
                    break

                assert line[0] == ">"
                transcript_id, chrom, region = line[1:].strip().split(":")
                sequence = transcriptome_file.readline().strip()
                transcripts[transcript_id] = sequence
        return transcripts

    def load_genome(self, file_path):
        """
        Read in a fasta file and load is a dictionary id -> sequence
        """
        genome = dict()
        with open(file_path) as transcriptome_file:
            while True:
                line = transcriptome_file.readline()
                if not line:
                    break

                assert line[0] == ">"
                contig = line[1:].strip()
                sequence = transcriptome_file.readline().strip()
                genome[contig] = sequence
        return genome

    def load_indels(self, file_path):
        """
        Read in the file of indel locations for a given custom genome

        Store it as a dictionary chrom -> (indel_starts, indel_data)
        where indel_starts is a numpy array of start locations of the indels
        (i.e. 1 based coordinates of the base in the custom genome where the insertion/deletion occurs immediately after)
        and indel_data is a list of tuples (indel_start, indel_type, indel_length)
        where indel_type is 'I' or 'D'

        The indel file is tab-separated with format "chrom:start type length"  and looks like the following:
        1:4897762       I       2
        1:7172141       I       2
        1:7172378       D       1

        Assumption is that the file is sorted by start and no indels overlap

        Moreover, return a dictionary chrom -> (offset_starts, offset_values)
        where offset_starts is as indel_starts, a sorted numpy array with indicating the positions
        where the offset (i.e. custom_genome_position - reference_genome_position values) change
        and offset_values is a list in the same order indicating these values.
        Note that the offset value is to be used for all bases AFTER the offset position, not on that base
        """
        indels = collections.defaultdict(lambda : (collections.deque(), collections.deque()))
        offset_values = collections.defaultdict(lambda : collections.deque())
        current_offsets = collections.defaultdict(lambda : 0)
        with open(file_path) as indel_file:
            for line in indel_file:
                loc, indel_type, length = line.split('\t')
                chrom, start = loc.split(':')

                # indel_start needs to be relative to the custom genome but the indel file has starts
                # that are relative to the reference genome
                indel_start = int(start) + current_offsets[chrom]
                indel_length = int(length)

                starts, data = indels[chrom]
                starts.append(indel_start)
                data.append((indel_start, indel_type, indel_length))

                if indel_type == "I":
                    current_offsets[chrom] += indel_length
                elif indel_type == "D":
                    current_offsets[chrom] -= indel_length

                offset_values[chrom].append(current_offsets[chrom])

        # We make these defaultdicts in case there are chromosomes without any indels
        # in which case we don't see them in this file, but we don't want to crash on them
        result = collections.defaultdict(lambda: (numpy.array([]), list()))
        offset_data = collections.defaultdict(lambda: (numpy.array([]), list()))
        for key, (starts, data) in indels.items():
            result[key] = (numpy.array(starts), list(data))
            offset_data[key] = (numpy.array(starts), offset_values[key])

        return result, offset_data

    def load_intron_quants(self, file_path):
        """
        Load an intron quantification file as two dictionaries,
        (transcript ID -> sum FPK of all introns in transcript) and
        (transcript ID -> list of FPKs of each intron in transcript)
        """
        transcript_intron_quants = dict() # Dictionary transcript -> FPK for all introns in the transcript, combined
        intron_quants = dict() # Dictioanry transcript -> array of FPKs for each intron in the transcript

        with open(file_path) as intron_quants_file:
            for line in intron_quants_file:

                if line.startswith("#"):
                    continue # Comment/header line

                gene, transcript, chrom, strand, transcript_intron_reads_FPK, intron_reads_FPK = line.strip().split("\t")

                transcript_intron_quants[transcript] = float(transcript_intron_reads_FPK)
                intron_quants[transcript] = [float(quant) for quant in intron_reads_FPK.split(",")]

        return transcript_intron_quants, intron_quants

    def load_gene_quants(self, file_path):
        """
        Read in a gene quantification file as two lists of gene IDs and of their read quantifications
        """
        genes = []
        gene_quants = []

        with open(file_path) as gene_quant_file:
            for line in gene_quant_file:
                if line.startswith("#"):
                    continue # Comment/header line

                gene, quant = line.strip().split("\t")

                genes.append(gene)
                gene_quants.append(float(quant))

        return genes, numpy.array(gene_quants)

    def load_isoform_quants(self, file_path):
        """
        Reads an isoform quant file into a dictionary gene -> (list of transcript IDs, list of psi values)
        """
        isoform_quants = dict()

        with open(file_path) as isoform_quant_file:
            for line in isoform_quant_file:

                if line.startswith("#"):
                    continue # Comment/header line

                gene, entries = line.strip().split("\t")
                isoforms = [entry.split(":") for entry in entries.split(",")]
                isoforms = [(isoform, float(psi)) for isoform, psi in isoforms]
                isoform_list = [isoform for isoform, psi in isoforms]
                psi_list = [psi for isoform, psi in isoforms]

                isoform_quants[gene] = (isoform_list, psi_list)

        return isoform_quants

    def load_allelic_quants(self, file_path):
        """
        Reads allelic quantification file into a dictionary: gene_id -> (allele 1 probability, allele 2 probability)
        """

        allelic_quant = dict()

        with open(file_path) as allele_quant_file:
            for line in allele_quant_file:
                if line.startswith("#"):
                    continue # Comment/header line

                gene, allele1, allele2 = line.split("\t")
                allele1 = float(allele1)
                allele2 = float(allele2)

                allelic_quant[gene] = (allele1, allele2)
        return allelic_quant

    def make_molecule(self):
        # Pick random gene
        gene_index = numpy.random.choice(len(self.genes), p=self.gene_probabilities)
        gene = self.genes[gene_index]
        gene_quant = self.gene_quants[gene_index]

        # Pick random transcript in gene
        transcripts, psis = self.isoform_quants[gene]
        transcript = numpy.random.choice(transcripts, p=psis)

        # Pick random allele based on the gene's allelic distribution
        allele_number = numpy.random.choice([1,2], p=self.allelic_quant[gene])

        # Read in annotation for the chosen transcript
        chrom,strand,tx_start,tx_end,starts,ends= self.annotations[allele_number - 1][transcript]

        # Determine if pre_mRNA or mature mRNA
        intron_quant = self.transcript_intron_quants[transcript]
        # TODO: check that this gives the appropriate fraction as pre_mRNA
        #       previously was using intron_quant / (intron_quant + gene_quant)
        #       but if assuming everything is either full pre_mRNA or mature mRNA then this should be
        #       the right fraction, which could happen to be greater than one (!)
        try:
            fraction_pre_mRNA = min(intron_quant / (gene_quant), 1)
        except ZeroDivisionError:
            # Should not ever get here since a gene with 0 gene_quant should have 0 chance of being chosen
            # however, if we do, we will just always give pre_mRNA
            fraction_pre_mRNA = 1.0

        pre_mRNA = numpy.random.uniform() < fraction_pre_mRNA
        if pre_mRNA:
            # If chosen to be pre_mRNA, overwrite the usual exon starts/ends with a single, big "exon"
            starts = [tx_start]
            ends = [tx_end]

        # Find cigar string relative to the reference genome (i.e. custom_genome_1 or custom_genome_2)
        gaps = [next_start - last_end - 1 for next_start,last_end in zip(starts[1:],ends[:-1])]
        cigar = ''.join( f"{end - start + 1}M{gap}N" for start,end,gap in zip(starts[:-1],ends[:-1],gaps)) \
                    + f"{ends[-1] - starts[-1] + 1}M"

        ref_cigar =  ''.join( f"{self.get_reference_cigar(start, end, chrom, allele_number)}{gap}N"
                            for start,end,gap in zip(starts[:-1],ends[:-1],gaps)) \
                        + self.get_reference_cigar(starts[-1], ends[-1], chrom, allele_number)
        ref_start = self.convert_genome_position_to_reference(starts[0], chrom, allele_number)

        transcript_id = f"{transcript}_{allele_number}{'_pre_mRNA' if pre_mRNA else ''}"


        # Build the actual sequence
        chrom_sequence = self.genomes[allele_number - 1][chrom]
        sequence = ''.join( chrom_sequence[start-1:end] for start,end in zip(starts, ends) )

        if strand == '-':
            # We always give the sequence from 5' to 3' end of the RNA molecule
            # so reverse complement this
            sequence = GeneralUtils.create_complement_strand(sequence)
            # NOTE: cigar string stays the same since that is relative to the + strand

        # TODO: for now, everything gets polyA but maybe shouldn't
        polyA_tail = True
        if polyA_tail:
            # TODO: polyA tails should vary in length
            # Add polyA tail to 3' end
            sequence = sequence + "A"*200
            # Soft-clip the polyA tail at the end since it shouldn't align
            if strand == "+":
                cigar = cigar + "200S"
                ref_cigar = ref_cigar + "200S"
            else:
                cigar = "200S" + cigar # Relative to + strand, the A's are going on the 5' end
                ref_cigar = "200S" + ref_cigar


        return sequence, starts[0], cigar, ref_start, ref_cigar, strand, chrom, transcript_id

    def get_reference_cigar(self, start, end, chrom, allele):
        """Returns the cigar string for the part of the custom chromosome on the
        segment from  start to end (inclusive, one based) relative to the reference
        genome
        """
        indel_starts, indel_data = self.indels[allele-1][chrom]

        cigar_components = []
        start_of_remaining = start

        # Skip to the first relevant indel, then continue until they go past this region
        index_idx = numpy.searchsorted(indel_starts, start)
        for indel in indel_data[index_idx:]:
            indel_start, indel_type, indel_length = indel
            if indel_start >= end:
                # The deletion or insertion begins immediately AFTER indel_start
                # so we're done, we've gone through all indels in our exon
                break

            if indel_start > start_of_remaining:
                # Match a region without indels
                cigar_components.append(f"{indel_start - start_of_remaining + 1}M")

            if indel_type == "I":
                # Custom genome has an insertion right after indel_start
                # Our region might end before the insertion does, so cap the length of the insertion
                length = min(end, indel_start + indel_length) - indel_start
                cigar_components.append(f"{length}I")
                # Increment to after the indel, which includes the inserted bases
                start_of_remaining = indel_start + indel_length + 1
            elif indel_type == "D":
                # Custom genome has a deletion right after indel_start
                cigar_components.append(f"{indel_length}D")
                # Advance to just after the deletion, but don't increment by the length of the deletion
                # since the deleted bases don't appear in our custom genome segment
                start_of_remaining = indel_start + 1

        # The last chunk of matches after the last indel
        if end >= start_of_remaining:
            cigar_components.append(f"{end - start_of_remaining + 1}M")

        return ''.join(cigar_components)

    def convert_genome_position_to_reference(self, position, chrom, allele):
        """
        Convert (1-indexed) position into the current (custom) genome into a position
        relative to the reference genome
        """
        # offset_positions is a sorted list of all offsets in this allele+chromosome
        # so find where our position fits in, so offset_index -1 is the last offset strictly before position
        offset_starts, offset_values = self.offset_data[allele-1][chrom]
        offset_index = numpy.searchsorted(offset_starts, position)

        if offset_index == 0:
            # If ours is before all offsets, then we need no offset
            return position
        else:
            offset = offset_values[offset_index-1]
            return position - offset

    def make_packet(self, sample, id="packet0", N=10_000):
        molecules = []
        for i in range(N):
            sequence, start, cigar, strand, ref_start, ref_cigar, chrom, transcript_id = self.make_molecule()
            mol = Molecule(
                    Molecule.new_id(transcript_id),
                    sequence,
                    start = start, # relative to the true ('parental') genome
                    cigar = cigar,
                    strand = strand,
                    source_start = ref_start, # relative to the reference genome
                    source_cigar = ref_cigar,
                    source_strand = strand,
                    transcript_id = transcript_id,
                    source_chrom = chrom)
            molecules.append(mol)
        return MoleculePacket(id, sample, molecules)

    def make_molecule_file(self, filepath, N=10_000):
        """
        Write out molecules to a tab-separated file

        Note: we write out a molecules start and cigar relative to the appropriate
        custom genome, either _1 or _2 as per the transcript id
        """
        with open(filepath, "w") as molecule_file:
            header = "#transcript_id\tchrom\tstart\tcigar\tref_start\tref_cigar\tstrand\tsequence\n"
            molecule_file.write(header)
            for i in range(N):
                sequence, start, cigar, ref_start, ref_cigar, strand, chrom, transcript_id = self.make_molecule()
                # NOTE: Not outputing the molecules start or cigar string since those are relative to parent
                #       which in this case is always trivial (start=1, cigar=###M) since the molecule is new
                line = "\t".join([transcript_id,
                                  chrom,
                                  str(start),
                                  cigar,
                                  str(ref_start),
                                  ref_cigar,
                                  strand,
                                  sequence]
                                  ) + "\n"

                molecule_file.write(line)

    def execute(self, sample, intron_quant_path, gene_quant_path, psi_quant_path,
                allele_quant_path, output_type, output_molecule_count, seed=None,
                molecules_per_packet=None):
        """This is the main method that generates simulated molecules and saves/
        exports them in the desired format. It uses the gene, transcript, intron,
        and allelic imbalance distributions generated by the other CAMPAREE steps.


        Parameters
        ----------
        sample : Sample
            Sample object corresponding to the input distributions. When exporting
            molecule packets, this Sample object is used to instantiate the
            MoleculePacket object.
        intron_quant_path : string
            Path to file containing intron expression distributions.
        gene_quant_path : string
            Path to file containing gene expression distributions.
        psi_quant_path : string
            Path to file containing transcript PSI value distributions.
        allele_quant_path : string
            Path to file containing distributions for allelic imbalance.
        output_type : string
            Type of file or object used to save or export simulated molecules.
            Sould be one of {', '.join(MoleculeMakerStep.OUTPUT_OPTIONS_W_EXTENSIONS.keys())}.
        output_molecule_count : integer
            Total number of molecules to save/export for the current Sample.
        seed : integer
            [OPTIONAL] Seed for random number generator. Used so repeated runs
            can produce the same results.
        molecules_per_packet : integer
            [OPTIONAL] Maximum number of molecules in each molecule packet. Must
            be positive, non-zero integer (this is not currently checked).

        """
        sample_data_directory = os.path.join(self.data_directory_path, f"sample{sample.sample_id}")
        log_file_path = os.path.join(self.log_directory_path, f'sample{sample.sample_id}',
                                     CAMPAREE_CONSTANTS.MOLECULE_MAKER_LOG_FILENAME)
        output_file_extension = MoleculeMakerStep.OUTPUT_OPTIONS_W_EXTENSIONS[output_type]

        if seed is not None:
            numpy.random.seed(seed)
        if not molecules_per_packet:
            molecules_per_packet=MoleculeMakerStep.DEFAULT_MOLECULES_PER_PACKET

        with open(log_file_path, "w") as log_file:

            print(f"Generating molecules for sample{sample.sample_id}.")
            log_file.write(f"Generating molecules for sample{sample.sample_id}.\n")

            log_file.write(f"Parameters:\n"
                           f"    Output file type: {output_type}\n"
                           f"    Output file extension: {output_file_extension}\n"
                           f"    Num molecules to generate: {output_molecule_count}\n"
                           f"    Num molecules per packet: {molecules_per_packet}\n"
                           f"    Random seed value: {seed}\n")

            log_file.write(f"Distribution files:\n"
                           f"    Intron quant file: {intron_quant_path}\n"
                           f"    Gene quant file: {gene_quant_path}\n"
                           f"    PSI quant file: {psi_quant_path}\n"
                           f"    Allele quant file: {allele_quant_path}\n")

            print('Loading gene, intron, transcript PSI, and allelic imbalance'
                  ' distributions.')
            log_file.write('Loading gene, intron, transcript PSI, and allelic'
                           ' imbalance distributions.\n')

            # Read and load data from gene, intron, transcript PSI, and allelic
            # imbalance distribution files.
            self.genes, self.gene_quants = self.load_gene_quants(gene_quant_path)
            self.gene_probabilities = self.gene_quants / numpy.sum(self.gene_quants)

            self.transcript_intron_quants, self.intron_quants = self.load_intron_quants(intron_quant_path)
            self.isoform_quants = self.load_isoform_quants(psi_quant_path)
            self.allelic_quant = self.load_allelic_quants(allele_quant_path)

            print('Loading annotations, transcriptome sequences, and genome sequences'
                  ' from botr parental genomes.')
            log_file.write('Loading annotations, transcriptome sequences, and genome'
                           ' sequences from both parental genomes.\n')

            # Read and load annotations, as well as full transcriptome and genome
            # sequences for each parental genome. This information is used when
            # generating the simulated molecule sequences.
            self.transcriptomes = \
                [self.load_transcriptome(os.path.join(sample_data_directory,
                                                      self._PARENTAL_TX_FASTA_FILENAME_PATTERN.format(genome_name=genome_name)))
                    for genome_name in [1,2]]
            self.annotations = \
                [self.load_annotation(os.path.join(sample_data_directory,
                                                   self._PARENTAL_ANNOT_FILENAME_PATTERN.format(genome_name=genome_name)))
                    for genome_name in [1,2]]
            self.genomes = \
                [self.load_genome(os.path.join(sample_data_directory,
                                               self._PARENTAL_GENOME_FASTA_FILENAME_PATTERN.format(genome_name=genome_name)))
                    for genome_name in [1,2]]

            print('Loading indel information from both parental genomes.')
            log_file.write('Loading indel information from both parental genomes.\n')

            # Read and load indel data for each parental genome. This information is
            # used when constructing CIGAR strings mapping transcripts back to their
            # locations in the original reference genome.
            indel_data = \
                [self.load_indels(os.path.join(sample_data_directory,
                                               self._PARENTAL_GENOME_INDEL_FILENAME_PATTERN.format(genome_name=genome_name)))
                    for genome_name in [1,2]]
            self.indels = [indels for indels, offset_data in indel_data]
            self.offset_data = [offset_data for indels, offset_data in indel_data]

            # Generate molecules and save/export them according to output type.
            print('Generating molecules and saving/exporting the results.')
            log_file.write('Generating molecules and saving/exporting the results.')
            if output_type == "packet":
                # TODO: potentially rounds down the number of molecules to make
                num_packets = output_molecule_count // molecules_per_packet
                for i in range(1,num_packets+1):
                    print(f"    Generating packet {i} of {num_packets}")
                    log_file.write(f"    Generating packet {i} of {num_packets}\n")
                    packet = self.make_packet(sample=sample, id=f"sample{sample.sample_id}.{i}")

                    molecule_packet_filename = os.path.join(sample_data_directory,
                                                            self.OUTPUT_FILENAME_PATTERN.format(output_type=output_type,
                                                                                                packet_num=i,
                                                                                                extension=output_file_extension))
                    with open(molecule_packet_filename, "wb") as out_file:
                        pickle.dump(packet, out_file)
            elif output_type == "file":
                print("Generating molecule file.")
                log_file.write("Generating molecule file.")
                molecule_output_filename = os.path.join(sample_data_directory,
                                                        self.OUTPUT_FILENAME_PATTERN.format(output_type=output_type,
                                                                                            packet_num="",
                                                                                            extension=output_file_extension))
                self.make_molecule_file(filepath=molecule_output_filename,
                                        N = output_molecule_count)

            log_file.write("\nALL DONE!\n")

    def get_commandline_call(self, sample, intron_quant_path, gene_quant_path,
                             psi_quant_path, allele_quant_path, output_type,
                             output_molecule_count, seed=None,
                             molecules_per_packet=None):
        """Prepare command to execute the MoleculeMakerStep from the command line,
        given all of the arugments used to run the execute() function.

        Parameters
        ----------
        sample : Sample
            Sample object corresponding to the input distributions. When exporting
            molecule packets, this Sample object is used to instantiate the
            MoleculePacket object.
        intron_quant_path : string
            Path to file containing intron expression distributions.
        gene_quant_path : string
            Path to file containing gene expression distributions.
        psi_quant_path : string
            Path to file containing transcript PSI value distributions.
        allele_quant_path : string
            Path to file containing distributions for allelic imbalance.
        output_type : string
            Type of file or object used to save or export simulated molecules.
            Sould be one of {', '.join(MoleculeMakerStep.OUTPUT_OPTIONS_W_EXTENSIONS.keys())}.
        output_molecule_count : integer
            Total number of molecules to save/export for the current Sample.
        seed : integer
            [OPTIONAL] Seed for random number generator. Used so repeated runs
            can produce the same results.
        molecules_per_packet : integer
            [OPTIONAL] Maximum number of molecules in each molecule packet. Must
            be positive, non-zero integer (this is not currently checked).

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.

        """
        #Retrieve path to the allelic_imbalance_quant.py script.
        molecule_maker_step_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        molecule_maker_step_path = molecule_maker_step_path.rstrip('c')

        command = (f" python {molecule_maker_step_path}"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --sample '{repr(sample)}'"
                   f" --intron_quant '{intron_quant_path}'"
                   f" --gene_quant '{gene_quant_path}'"
                   f" --psi_quant '{psi_quant_path}'"
                   f" --allele_quant '{allele_quant_path}'"
                   f" --output_type {output_type}"
                   f" --output_molecule_count {output_molecule_count}")
        if seed is not None:
            command += f" --seed {seed}"
        if molecules_per_packet:
            command += f" --molecules_per_packet {molecules_per_packet}"

        return command

    def get_validation_attributes(self, sample, intron_quant_path, gene_quant_path,
                                  psi_quant_path, allele_quant_path, output_type,
                                  output_molecule_count, seed=None,
                                  molecules_per_packet=None):
        """Prepare attributes required by is_output_valid() function to validate
        output generated by the MoleculeMakerStep job.

        Parameters
        ----------
        sample : Sample
            Sample object corresponding to the input distributions. When exporting
            molecule packets, this Sample object is used to instantiate the
            MoleculePacket object.
        intron_quant_path : string
            Path to file containing intron expression distributions. [Note: this
            parameter is captured just so get_validation_attributes() accepts
            the same arguments as get_commandline_call(). It is not used here.]
        gene_quant_path : string
            Path to file containing gene expression distributions. [Note: this
            parameter is captured just so get_validation_attributes() accepts
            the same arguments as get_commandline_call(). It is not used here.]
        psi_quant_path : string
            Path to file containing transcript PSI value distributions. [Note:
            this parameter is captured just so get_validation_attributes()
            accepts the same arguments as get_commandline_call(). It is not used
            here.]
        allele_quant_path : string
            Path to file containing distributions for allelic imbalance. [Note:
            this parameter is captured just so get_validation_attributes()
            accepts the same arguments as get_commandline_call(). It is not used
            here.]
        output_type : string
            Type of file or object used to save or export simulated molecules.
            Sould be one of {', '.join(MoleculeMakerStep.OUTPUT_OPTIONS_W_EXTENSIONS.keys())}.
        output_molecule_count : integer
            Total number of molecules to save/export for the current Sample.
        seed : integer
            [OPTIONAL] Seed for random number generator. Used so repeated runs
            can produce the same results. [Note: this parameter is captured just
            so get_validation_attributes() accepts the same arguments as
            get_commandline_call(). It is not used here.]
        molecules_per_packet : integer
            [OPTIONAL] Maximum number of molecules in each molecule packet. Must
            be positive, non-zero integer (this is not currently checked).

        Returns
        -------
        dict
            A MoleculeMakerStep job's data_directory, log_directory, corresponding
            sample ID, output file type, output molecule count, and the number of
            molecules per packet.
        """
        validation_attributes = {}
        validation_attributes['data_directory'] = self.data_directory_path
        validation_attributes['log_directory'] = self.log_directory_path
        validation_attributes['sample_id'] = sample.sample_id
        validation_attributes['output_type'] = output_type
        validation_attributes['output_molecule_count'] = output_molecule_count
        validation_attributes['molecules_per_packet'] = molecules_per_packet
        return validation_attributes

    @staticmethod
    def is_output_valid(validation_attributes):
        """Check if output of MoleculeMakerStep for a specific job/execution is
        correctly formed and valid, given a job's data directory, log directory,
        sample ID, output file type, output molecule count, and the number of
        molecules per packet (if provided). Prepare these attributes for a given
        job using the get_validation_attributes() method.

        Parameters
        ----------
        validation_attributes : dict
            A job's data_directory, log_directory, corresponding sample_id,
            output file type, output molecule count, and the number of molecules
            per packet (if provided).

        Returns
        -------
        boolean
            True  - MoleculeMakerStep output files were created and are well formed.
            False - MoleculeMakerStep output files do not exist or are missing data.

        """

        data_directory_path = validation_attributes['data_directory']
        log_directory_path = validation_attributes['log_directory']
        sample_id = validation_attributes['sample_id']
        output_type = validation_attributes['output_type']
        output_molecule_count = validation_attributes['output_molecule_count']
        molecules_per_packet = validation_attributes.get('molecules_per_packet')
        output_file_extension = MoleculeMakerStep.OUTPUT_OPTIONS_W_EXTENSIONS[output_type]

        if not molecules_per_packet:
            molecules_per_packet = MoleculeMakerStep.DEFAULT_MOLECULES_PER_PACKET

        valid_output = False

        # Construct output filenames/paths
        sample_data_directory = os.path.join(data_directory_path, f"sample{sample_id}")
        log_file_path = os.path.join(log_directory_path, f'sample{sample_id}',
                                     CAMPAREE_CONSTANTS.MOLECULE_MAKER_LOG_FILENAME)

        # Check existence of output files. If output type is molecule packet,
        # there will be multiple output files and this code will check that they
        # all exist (based on number of molecules per packet).
        all_molecule_files_exist = False
        if output_type == "packet":
            total_num_packets = output_molecule_count // molecules_per_packet
            num_packets_exist = 0
            for i in range(1, total_num_packets+1):
                molecule_packet_file = os.path.join(sample_data_directory,
                                                    MoleculeMakerStep.OUTPUT_FILENAME_PATTERN.format(output_type=output_type,
                                                                                                     packet_num=i,
                                                                                                     extension=output_file_extension))
                if os.path.isfile(molecule_packet_file):
                    num_packets_exist += 1

            if total_num_packets == num_packets_exist:
                all_molecule_files_exist = True

        elif output_type == "file":
            molecule_file = os.path.join(sample_data_directory,
                                         MoleculeMakerStep.OUTPUT_FILENAME_PATTERN.format(output_type=output_type,
                                                                                          packet_num="",
                                                                                          extension=output_file_extension))
            all_molecule_files_exist = os.path.isfile(molecule_file)


        # TODO: Report reason why out validation failed.

        if all_molecule_files_exist and os.path.isfile(log_file_path):

            #Read last line in log file
            line = ""
            with open(log_file_path, "r") as log_file:
                for line in log_file:
                    line = line.rstrip()
            if line == "ALL DONE!":
                valid_output = True

        return valid_output

    @staticmethod
    def main():
        """Entry point into script. Parses the argument list to obtain all the
        files needed and feeds them to the class constructor. Calls the appropriate
        methods thereafter.

        """
        parser = argparse.ArgumentParser(description='Generate simulated molecules and export/save.')
        parser.add_argument('-l', '--log_directory_path', required=True,
                            help="Path to log directory.")
        parser.add_argument('-d', '--data_directory_path', required=True,
                            help='Path to data directory')
        parser.add_argument('--sample', required=True,
                            help='String representation of a Sample object.')
        parser.add_argument('--intron_quant', required=True,
                            help='Path to file containing intron expression distributions.')
        parser.add_argument('--gene_quant', required=True,
                            help='Path to file containing gene expression distributions.')
        parser.add_argument('--psi_quant', required=True,
                            help='Path to file containing transcript PSI value distributions.')
        parser.add_argument('--allele_quant', required=True,
                            help='Path to file containing distributions for allelic imbalance.')
        parser.add_argument('--output_type', required=True,
                            help=f"Type of molecule output ({', '.join(MoleculeMakerStep.OUTPUT_OPTIONS_W_EXTENSIONS.keys())}).")
        parser.add_argument('--output_molecule_count', type=int, required=True,
                            help='Number of molecules to generate.')
        parser.add_argument('--seed', type=int, default=None, required=False,
                            help='Seed value for random number generator.')
        parser.add_argument('--molecules_per_packet', type=int, default=None, required=False,
                            help='Number of molecules per molecule packet. '
                                 'Only used if output_type set to "packet".')
        args = parser.parse_args()
        sample = eval(args.sample)

        molecule_maker = MoleculeMakerStep(log_directory_path=args.log_directory_path,
                                           data_directory_path=args.data_directory_path)
        molecule_maker.execute(sample=sample,
                               intron_quant_path=args.intron_quant,
                               gene_quant_path=args.gene_quant,
                               psi_quant_path=args.psi_quant,
                               allele_quant_path=args.allele_quant,
                               output_type=args.output_type,
                               output_molecule_count=args.output_molecule_count,
                               seed=args.seed,
                               molecules_per_packet=args.molecules_per_packet)

if __name__ == "__main__":
    sys.exit(MoleculeMakerStep.main())
