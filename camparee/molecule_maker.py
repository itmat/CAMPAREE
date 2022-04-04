import os
import json
import pathlib
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
from beers_utils.cigar import chain_from_splits, split_cigar
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

    def __init__(self, log_directory_path, data_directory_path=None, parameters=None):
        """Constructor for MoleculeMakerStep object.

        Parameters
        ----------
        log_directory_path : string
            Full path to log directory.
        parameters : dict
            Dictionary of other parameters specified by the config file. This
            parameter is not used by this class and is retained for uniformity
            with all other CAMPAREE steps.

        """
        self.log_directory_path = log_directory_path
        self.data_directory_path = data_directory_path
        self.min_polyA_tail_length = parameters["min_polyA_tail_length"]
        self.max_polyA_tail_length = parameters["max_polyA_tail_length"]
        self.parameters = parameters

    # Nearly all of the validation for this step is already performed in the
    # expression_pipeline, since the output options in the config file are
    # specified outside of the standard 'steps' framework.
    def validate(self):
        if (self.min_polyA_tail_length < 0):
            return False
        if (self.min_polyA_tail_length > self.max_polyA_tail_length):
            return False
        if (self.max_polyA_tail_length < 0):
            return False
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

        NOTE: adds a 'padding' M (match) of length 1_000_000_000 to the end of the genome
        since we don't load the full (reference) genome length in but may still need to
        get locations past the last indel
        """
        genome_cigars = collections.defaultdict(lambda : collections.deque())
        last_indexes = collections.defaultdict(lambda : 1)
        with open(file_path) as indel_file:
            for line in indel_file:
                loc, indel_type, length = line.split('\t')
                chrom, start = loc.split(':')

                start = int(start)
                length = int(length)

                last_index = last_indexes[chrom]
                if last_index < start:
                    # Match up to the start of the indel
                    genome_cigars[chrom].append(
                        ('M', start - last_index + 1)
                    )

                # Add the indel
                if indel_type == 'I':
                    genome_cigars[chrom].append(
                        ('I', length)
                    )
                else:
                    genome_cigars[chrom].append(
                        ('D', length)
                    )

        # Gather into a results
        # which is a default dict which fills in the padding to
        # any chromosome even if it doesn't have any indels
        PADDING = ("M", 1_000_000_000)
        results = collections.defaultdict(lambda : [PADDING])
        for chrom, split_cigar in genome_cigars.items():
            results[chrom] = list(split_cigar)

        return results

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

    def make_molecule(self, sample):
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

        cigar_split = split_cigar(cigar)
        ref_start, ref_cigar, _ = chain_from_splits(
                starts[0], cigar_split, strand,
                1, self.genome_cigar_splits[allele_number - 1][chrom], "+"
        )

        transcript_id = f"{sample.sample_id}_{transcript}_{allele_number}{'_pre_mRNA' if pre_mRNA else ''}"


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
            polyA_length = numpy.random.randint(self.min_polyA_tail_length, self.max_polyA_tail_length + 1)
            sequence = sequence + "A"*polyA_length
            # Soft-clip the polyA tail at the end since it shouldn't align
            if strand == "+":
                cigar = cigar + f"{polyA_length}S"
                ref_cigar = ref_cigar + f"{polyA_length}S"
            else:
                cigar =   f"{polyA_length}S" + cigar # Relative to + strand, the A's are going on the 5' end
                ref_cigar = f"{polyA_length}S" + ref_cigar


        return sequence, starts[0], cigar, ref_start, ref_cigar, strand, chrom, transcript_id

    def make_packet(self, sample, id="packet0", N=10_000):
        molecules = []
        for i in range(N):
            sequence, start, cigar, ref_start, ref_cigar, strand, chrom, transcript_id = self.make_molecule(sample)
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

    def make_molecule_file(self, filepath, sample, N=10_000):
        """
        Write out molecules to a tab-separated file

        Note: we write out a molecules start and cigar relative to the appropriate
        custom genome, either _1 or _2 as per the transcript id
        """
        with open(filepath, "w") as molecule_file:
            header = "#transcript_id\tchrom\tstart\tcigar\tref_start\tref_cigar\tstrand\tsequence\n"
            molecule_file.write(header)
            for i in range(N):
                sequence, start, cigar, ref_start, ref_cigar, strand, chrom, transcript_id = self.make_molecule(sample)
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

    def execute(self, sample, sample_data_directory, output_type, output_molecule_count, seed=None,
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
        sample_data_directory : string
            Path to directory containing the data for the sample.
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
        sample_log_dir = pathlib.Path(self.log_directory_path) / f'sample{sample.sample_id}'
        sample_log_dir.mkdir(exist_ok=True)
        log_file_path = sample_log_dir / CAMPAREE_CONSTANTS.MOLECULE_MAKER_LOG_FILENAME
        output_file_extension = MoleculeMakerStep.OUTPUT_OPTIONS_W_EXTENSIONS[output_type]

        if seed is not None:
            numpy.random.seed(seed)
        if not molecules_per_packet:
            molecules_per_packet=MoleculeMakerStep.DEFAULT_MOLECULES_PER_PACKET

        with open(log_file_path, "w") as log_file:

            print(f"Generating molecules for sample{sample.sample_id}.")
            log_file.write(f"Generating molecules for sample{sample.sample_id}.\n")

            intron_quant_path = os.path.join(sample_data_directory, CAMPAREE_CONSTANTS.INTRON_OUTPUT_FILENAME)
            gene_quant_path = os.path.join(sample_data_directory, CAMPAREE_CONSTANTS.TXQUANT_OUTPUT_GENE_FILENAME)
            psi_quant_path = os.path.join(sample_data_directory, CAMPAREE_CONSTANTS.TXQUANT_OUTPUT_PSI_FILENAME)
            allele_quant_path = os.path.join(sample_data_directory, CAMPAREE_CONSTANTS.ALLELIC_IMBALANCE_OUTPUT_FILENAME)

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
                  ' from both parental genomes.')
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
            self.genome_cigar_splits =  [self.load_indels(os.path.join(sample_data_directory,
                                               self._PARENTAL_GENOME_INDEL_FILENAME_PATTERN.format(genome_name=genome_name)))
                                            for genome_name in [1,2]]

            # Generate molecules and save/export them according to output type.
            print('Generating molecules and saving/exporting the results.')
            log_file.write('Generating molecules and saving/exporting the results.')
            print(f"Molecule maker output type {repr(output_type)}")
            if output_type == "packet":
                # TODO: potentially rounds down the number of molecules to make
                num_packets = output_molecule_count // molecules_per_packet
                for i in range(1,num_packets+1):
                    print(f"    Generating packet {i} of {num_packets}")
                    log_file.write(f"    Generating packet {i} of {num_packets}\n")
                    packet = self.make_packet(sample=sample, id=f"sample{sample.sample_id}.{i}", N=molecules_per_packet) #TODO: id needs to be an integer

                    molecule_packet_filename = os.path.join(sample_data_directory,
                                                            self.OUTPUT_FILENAME_PATTERN.format(output_type=output_type,
                                                                                                packet_num=i,
                                                                                                extension=output_file_extension))
                    with open(molecule_packet_filename, "wb") as out_file:
                        pickle.dump(packet, out_file)
            elif output_type == "file":
                molecule_output_filename = os.path.join(sample_data_directory,
                                                        self.OUTPUT_FILENAME_PATTERN.format(output_type=output_type,
                                                                                            packet_num="",
                                                                                            extension=output_file_extension))
                print(f"Generating molecule file {molecule_output_filename}.")
                log_file.write(f"Generating molecule file {molecule_output_filename}.")
                self.make_molecule_file(filepath=molecule_output_filename,
                                        N = output_molecule_count,
                                        sample = sample)
            elif output_type == "generator":
                def generator():
                    num_packets = output_molecule_count // molecules_per_packet
                    print(f"Generating {num_packets} packets")
                    for i in range(1, num_packets+1):
                        packet = self.make_packet(sample=sample, id=i, N=molecules_per_packet)
                        yield packet
                return generator()
            else:
                raise ValueError(f"Expected output_type to be 'packet', 'file', or 'generator'. Instead got {repr(output_type)}")

            log_file.write("\nALL DONE!\n")

    def get_commandline_call(self, sample, sample_data_directory,
                             output_type, output_molecule_count,
                             seed=None,
                             molecules_per_packet=None):
        """Prepare command to execute the MoleculeMakerStep from the command line,
        given all of the arugments used to run the execute() function.

        Parameters
        ----------
        sample : Sample
            Sample object corresponding to the input distributions. When exporting
            molecule packets, this Sample object is used to instantiate the
            MoleculePacket object.
        sample_data_directory : string
            Path to directory containing the sample data
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
                   f" --parameters '{json.dumps(self.parameters)}'"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --sample_data_directory {sample_data_directory}"
                   f" --sample '{repr(sample)}'"
                   f" --output_type {output_type}"
                   f" --output_molecule_count {output_molecule_count}")
        if seed is not None:
            command += f" --seed {seed}"
        if molecules_per_packet:
            command += f" --molecules_per_packet {molecules_per_packet}"

        return command

    def get_validation_attributes(self, sample, sample_data_directory,
                                  output_type,
                                  output_molecule_count,
                                  seed=None,
                                  molecules_per_packet=None):
        """Prepare attributes required by is_output_valid() function to validate
        output generated by the MoleculeMakerStep job.

        Parameters
        ----------
        sample : Sample
            Sample object corresponding to the input distributions. When exporting
            molecule packets, this Sample object is used to instantiate the
            MoleculePacket object.
        sample_data_path : string
            Path to directory containing all the sample data.
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
            A MoleculeMakerStep job's sample_data_directory, log_directory, corresponding
            sample ID, output file type, output molecule count, and the number of
            molecules per packet.
        """
        validation_attributes = {}
        validation_attributes['log_directory'] = self.log_directory_path
        validation_attributes['sample_data_directory'] = sample_data_directory
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

        sample_data_directory = validation_attributes['sample_data_directory']
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
        elif output_type == "generator":
            pass # Always valid
        else:
            # Unknown output type specified
            valid_output = False


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
        parser.add_argument('--parameters', required=True,
                            help="JSON of parameters")
        parser.add_argument('-d', '--sample_data_directory', required=True,
                            help='Path to sample data directory')
        parser.add_argument('--sample', required=True,
                            help='String representation of a Sample object.')
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

        molecule_maker = MoleculeMakerStep(
                log_directory_path=args.log_directory_path,
                parameters = json.loads(args.parameters))
        molecule_maker.execute(sample=sample,
                               sample_data_directory=args.sample_data_directory,
                               output_type=args.output_type,
                               output_molecule_count=args.output_molecule_count,
                               seed=args.seed,
                               molecules_per_packet=args.molecules_per_packet)

if __name__ == "__main__":
    sys.exit(MoleculeMakerStep.main())
