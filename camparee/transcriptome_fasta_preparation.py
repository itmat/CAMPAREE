import argparse
import re
import sys
import os

from camparee.abstract_camparee_step import AbstractCampareeStep
from camparee.camparee_constants import CAMPAREE_CONSTANTS

class TranscriptomeFastaPreparationStep(AbstractCampareeStep):
    """Produces a transcriptome FASTA file, given a genome FASTA file, a file
    containing exon locations, and an annotation file. Additionally, any line in
    the annotation file related to a chromosome not available in the genome fasta
    file is discarded in a new, trimmed version of the annotation file.

    The object is constructed with 2 input file sources (genome fasta, annotation)
    and 2 output file sources (trimmed annotation, transcriptome fasta).
    Additionally another output file, named like the genome fasta file but
    suffixed with '_edited' contains a munged version of the genome fasta file
    where each chromosome sequence occupies one line.

    """

    Name = "Transcriptome FASTA Preparation Step"

    def __init__(self, log_directory_path, data_directory_path, parameters = dict()):
        """Constructor for TranscriptomeFastaPreparationStep object.

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
        # TODO: Could probably change to a class variable, or a constant contained in another class.
        # Provides the regex pattern for extracting components of the exon location string
        self.exon_info_pattern = re.compile(r'(.*):(\d+)-(\d+)')


    def validate(self):
        return True

    def execute(self, sample_id, genome_suffix, genome_fasta_file_path, annotation_file_path,
                include_suffix_w_tx_id=False):
        """Main work-horse function that does the work of creating a transcriptome
        fasta file from the provided inputs.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to this reference genome. Used to
            construct output and log paths for this specific execution.
        genome_suffix : string
            Suffix to identify the parent/allele of the source genome. Should be
            1 or 2. This same suffix is a appended to all output files, and the
            individual transcript IDs in the transcriptome FASTA if the
            include_suffix_w_tx_id parameter is set to TRUE.
        genome_fasta_file_path : string
            Input genome fasta filename containing all the chromosomes of interest.
            No line breaks are allowed within the chromosome sequence. This is
            generally prepared by the GenomeBuilderStep, in which case it should
            have no line breaks within the chromosome sequence.
        annotation_file_path : string
            Input transcript annotation file - fields are (chromosome, strand,
            start, end, exon count, exon starts, exon ends, transcript ID, etc.).
            This is generally prepared by the UpdateAnnotationForGenomeStep.
        include_suffix_w_tx_id : boolean
            Append parent/allele suffix to transcript names in FASTA headers of
            the output file when set to True [Default: False].

        """

        self.sample_id = sample_id
        self.genome_suffix = genome_suffix
        self.genome_fasta_file_path = genome_fasta_file_path
        self.edited_genome_fasta_file_path = os.path.splitext(genome_fasta_file_path)[0] + "_edited.fa"
        self.annotation_file_path = annotation_file_path
        self.trimmed_annotation_file_path = os.path.splitext(annotation_file_path)[0] + "_trimmed.txt"
        self.include_suffix_w_tx_id = include_suffix_w_tx_id

        self.transcriptome_fasta_file_path = os.path.join(self.data_directory_path, f'sample{self.sample_id}',
                                                          CAMPAREE_CONSTANTS.TRANSCRIPTOME_FASTA_OUTPUT_FILENAME_PATTERN.format(genome_name=genome_suffix))
        self.log_file_path = os.path.join(self.log_directory_path, f'sample{self.sample_id}',
                                          CAMPAREE_CONSTANTS.TRANSCRIPTOME_FASTA_LOG_FILENAME_PATTERN.format(genome_name=genome_suffix))

        # Holds unique listing of exon locations
        self.exon_location_list = set()
        # Dictionaries to record whether a chromosome is available in the genome file, is available in the exon
        # file.
        self.chromosome_in_genome_file = dict()
        self.chromosome_in_exon_file = dict()
        # Since the genes fasta file will be open for appending, we need to insure that the file doesn't
        # currently exist.
        try:
            os.remove(self.transcriptome_fasta_file_path)
        except OSError:
            pass

        with open(self.log_file_path, 'w') as log_file:

            # THIS IS NOT CURRENTLY NECESSARY, AS THE PRECEDING STEPS CREATE A GENOME
            # FASTA IN THE APPROPRIATE FORMAT. ALTERNATIVELY, A FUNCTION THAT PERFORMS
            # THESE OPERATIONS IS AVAILABLE IN CampareeUtils.edit_reference_genome().
            # NOTE: THE FUNCTION HERE PROCESSES THE FASTA FILE AS IT'S READ AND DOES
            #       NOT NEED TO STORE THE ENTIRE FILE IN MEMORY. WE SHOULD UPDATE THE
            #       FUNCTION IN CampareeUtils, TO OPERATE LIKE THIS AS WELL.
            # Munge the original genome fasta file to create an edited version in which the chromosome sequence has
            # no internal line breaks and where all bases are in upper case.
            self.scrub_genome_fasta_file()
            log_file.write(f"done scrubbing genome fasta file\n")

            # Create a unique listing of exon locations from the annotation file data.
            self.create_exon_location_list()
            log_file.write(f"done identifying unique set of exon locations\n")

            log_file.write(f"Assembling sequences for each transcript by chromosome\n")

            # Open the genome fasta file for reading only.
            with open(self.edited_genome_fasta_file_path, 'r') as genome_fasta_file:

                # Iterate over each chromosome in the genome fasta file.  Note that
                # the chromosome sequence is expected on only one line at this stage.
                for line in genome_fasta_file:

                    # Remove the leading '>' character to leave the chromosome
                    chromosome = line.lstrip('>').rstrip('\n')

                    # Note that the chromosome 'chromosome' is among those listed in
                    # the genome file
                    self.chromosome_in_genome_file[chromosome] = True

                    # Collect the chromosome sequence from the following line
                    sequence = genome_fasta_file.readline().rstrip('\n')

                    # Use the chromosome and its sequence to construct a mapping of
                    # exon location string to their sequences.
                    exon_sequence_map = self.create_exon_sequence_map(chromosome, sequence)
                    log_file.write(f"\tdone with exons for {chromosome}\n")

                    # Generate the transcript fasta file using the genome chromosome
                    # and the related exon sequence map.
                    self.make_tx_fasta_file(chromosome, exon_sequence_map)
                    log_file.write(f"\tdone with transcripts for {chromosome}\n")

        # Finally create a trimmed annotation file with any transcripts unrelated
        # to the given genome chromosomes provided, discarded.
        # This method saves messages to the log file internally, so it's outside
        # of the with block above to avoid opening the log file from multiple
        # sources.
        self.trim_annotation_file()

        #Final entry in log file to indicate execution() method finished.
        with open(self.log_file_path, 'a') as log_file:
            log_file.write('\nALL DONE!\n')

    def scrub_genome_fasta_file(self):
        """
        Edits the genome fasta file, creating an edited version (genome fasta filename without extension + _edited.fa).
        Edits include:
        1.  Removing suplemmental information from the description line
        2.  Removing internal newlines in the sequence
        3.  Insuring all bases in sequence are represented in upper case.
        This edited file is the one used in subsequent scripts.
        """

        # A flag to denote when a fasta sequence is being processed
        in_sequence = False

        # Open the original genome fasta file for reading only and a new edited version of the genome fasta
        # file for writing only.
        with open(self.genome_fasta_file_path, 'r') as genome_fasta_file, \
                open(self.edited_genome_fasta_file_path, 'w') as edited_genome_fasta_file:

            # Iterate over the lines in the original genome fasta file
            for line in genome_fasta_file:

                # Identify whether the current line is a description line or a sequence line
                if line.startswith('>'):

                    # For a description line, remove any supplemental information following the identifier and
                    # if the in_sequence flag is raised, lower it and add a line break to the new genome fasta
                    # file before adding the modified description line.
                    identifier_only = re.sub(r'[ \t].*', '', line)
                    if in_sequence:
                        edited_genome_fasta_file.write("\n")
                        in_sequence = False
                    edited_genome_fasta_file.write(identifier_only)
                # Otherwise, add the sequence to the new genome fasta file after removing the line break and
                # insuring all bases are in upper case.  Also raise the in sequence flag.
                else:
                    edited_genome_fasta_file.write(line.rstrip('\n').upper())
                    in_sequence = True

            # Finally add a line break to the end of the new genome fasta file.
            edited_genome_fasta_file.write("\n")

    def create_exon_location_list(self):
        """
        Generate a unique listing of exon location strings from the provided
        annotation file.  Note that the same exon may appear in multiple transcripts.
        So the listing is actually a set to avoid duplicate entries.
        """

        # Open the annotation file for reading only
        with open(self.annotation_file_path, 'r') as annotation_file:

            # Iterate over each line in the file
            for line in annotation_file:
                if line.startswith("#"): #Comment line
                    continue

                # Collect all the field values for the line read following newline removal.
                (chromosome, strand, start, end, exon_count, exon_starts, exon_ends, name, *other) = \
                    line.rstrip('\n').split('\t')

                # Remove any trailing commas in the exon starts and exon ends fields
                # and split the starts and stops into their corresponding lists
                exon_starts_list = re.sub(r'\s*,\s*$', '', exon_starts).split(",")
                exon_ends_list = re.sub(r'\s*,\s*$', '', exon_ends).split(",")

                # For each exon belonging to the transcript, construct the exon's
                # location string and add it to the growing set of exon locations
                # assuming it is not already present.
                for index in range(int(exon_count)):
                    exon_location = f'{chromosome}:{(int(exon_starts_list[index]) + 1)}-{(exon_ends_list[index])}'
                    self.exon_location_list.add(exon_location)

    def trim_annotation_file(self):
        """
        Create a trimmed annotation file in which lines related to chromosomes
        that are not present in the genome fasta file are expunged.  If there
        are no such omissions, the files will be identical.
        """

        # TODO: rather than always create this file, it might be better if it's
        #       only created when missing_genome_chromosomes is non-empty. Then
        #       any future steps that come looking for this file will check to
        #       see if the edited one exists, and load the original annotation if
        #       it does not.

        # Holds a list of chromosomes found for exons in the annotation file that are not available in the
        # genome fasta file
        missing_genome_chromosomes = []

        with open(self.log_file_path, 'a') as log_file:

            log_file.write(f"Identifying chromosomes in the annotation with no corresponding genome sequence.\n")

            # Iterate over all the chromosomes found for exons in the annotation file
            for chromosome in self.chromosome_in_exon_file.keys():

                # If that chromosome is not also in the genome fasta file, issue a warning and add the
                # unavailable chromosome to the list of missing chromosomes.  Otherwise, note the
                # chromosome as available.
                if chromosome not in self.chromosome_in_genome_file:
                    log_file.write(f"\tno genome sequence for {chromosome}\n")
                    missing_genome_chromosomes.append(chromosome)
                else:
                    log_file.write(f"\tsequence available for {chromosome}\n")

            # If there are missing chromosomes, note that fact.
            if missing_genome_chromosomes:
                log_file.write(f"Removing the transcripts on chromosomes for which no genome sequences are available.\n")

        # Open the original annotation for reading and the edited annotation file for writing.
        with open(self.annotation_file_path, 'r') as annotation_in, \
                open(self.trimmed_annotation_file_path, 'w') as annotation_out:

            # Iterate over the lines in the original annotation file
            for line in annotation_in:

                # If none of the missing chromosomes are found on the line, write
                # the line out to the trimmed version of the annotation file.
                if not any(chromosome in line for chromosome in missing_genome_chromosomes):
                    annotation_out.write(line)

    def create_exon_sequence_map(self, genome_chromosome, sequence):
        """
        For the given genome chromosome and its sequence, create a dictionary of exon sequences keyed to the exon's
        location (i.e., chr:start-end).
        :param genome_chromosome: given genome chromosome
        :param sequence: the genome sequence corresponding to the genome chromosome (without line breaks)
        :return: map of exon location : exon sequence
        """

        # Start with a empty dictionary
        exon_sequence_map = dict()

        # Iterate over the list of previously obtained exon locations
        for exon_location in self.exon_location_list:

            # Extract the chromosome, start and end from the exon location string
            exon_info_match = re.search(self.exon_info_pattern, exon_location)
            chromosome = exon_info_match.group(1)
            exon_start = int(exon_info_match.group(2))
            exon_end = int(exon_info_match.group(3))

            # Note that the exon listing contains at least one exon on chromosome 'chromosome'
            self.chromosome_in_exon_file[chromosome] = True

            # If the chromosome on which the exon is located is the same of the genome chromosome
            # provided as a parameter, get the sequence for that exon and create a dictionary
            # entry relating the exon location string to the exon sequence.
            if chromosome == genome_chromosome:
                exon_sequence_map[exon_location] = sequence[exon_start-1:exon_end]

        # Return the dictionary of exon location string : exon sequence for the genome chromosome
        # provided in the parameter list.
        return exon_sequence_map

    def make_tx_fasta_file(self, genome_chromosome, exon_sequence_map):

        # Open the original annotation file for reading and the genes fasta file for appending.
        with open(self.annotation_file_path, 'r') as annotation_file, \
             open(self.transcriptome_fasta_file_path, 'a') as transcriptome_fasta_file:

            # Iterate over the annotation file
            for line in annotation_file:
                if line.startswith("#"): # Comment line
                    continue

                # Collect all the field values for the line read following newline removal.
                (chromosome, strand, start, end, exon_count, exon_starts, exon_ends, feature_id, *other) =\
                    line.rstrip('\n').split('\t')

                # Remove trailing commas in the exon starts and exon ends fields and split
                # the starts and stops into their corresponding lists
                exon_starts_list = re.sub(r'\s*,\s*$', '', exon_starts).split(",")
                exon_ends_list = re.sub(r'\s*,\s*$', '', exon_ends).split(",")

                # If the chromosome on which this transcript is located is the
                # same of the genome chromosome provided as a parameter, render
                # this transcript as an entry in the genes fasta file
                if chromosome == genome_chromosome:

                    # Initialize the transcripts's sequence
                    tx_sequence = ""

                    # For each exon belonging to the transcript, construct the
                    # exon's location string and use it as a key to obtain the
                    # actual exon sequence.  Concatenate that exon sequence to
                    # the transcript sequence.  Note that the 1 added to each
                    # exon start takes into account the zero based and half-
                    # open ucsc coordinates.
                    for index in range(int(exon_count)):
                        exon_key = f'{chromosome}:{(int(exon_starts_list[index]) + 1)}-{(exon_ends_list[index])}'
                        exon_sequence = exon_sequence_map[exon_key]
                        tx_sequence += exon_sequence

                    # Remove some spurious characters from the transcript ID
                    tx_id = re.sub(r'::::.*', '', feature_id)
                    # TODO determine if this substitution is still needed.
                    tx_id = re.sub(r'\([^(]+$', '', tx_id)

                    if self.include_suffix_w_tx_id:
                        tx_id = tx_id + "_" + self.genome_suffix

                    # Write the 1st line of the fasta entry - transcript location string
                    transcriptome_fasta_file.write(f'>{tx_id}:{chromosome}:{start}-{end}_{strand}\n')

                    # Write the 2nd line of the fasta entry - gene sequence.
                    transcriptome_fasta_file.write(tx_sequence + '\n')

    def get_commandline_call(self, sample_id, genome_suffix, genome_fasta_file_path,
                             annotation_file_path, include_suffix_w_tx_id=False):
        """
        Prepare command to execute the TranscriptomeFastaPreparationStep from
        the command line, given all of the arugments used to run the execute()
        function.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to this reference genome. Used to
            construct output and log paths for this specific execution.
        genome_suffix : string
            Suffix to identify the parent/allele of the source genome. Should be
            1 or 2. This same suffix is a appended to all output files, and the
            individual transcript IDs in the transcriptome FASTA if the
            include_suffix_w_tx_id parameter is set to TRUE.
        genome_fasta_file_path : string
            Input genome fasta filename containing all the chromosomes of interest.
            No line breaks are allowed within the chromosome sequence.
        annotation_file_path : string
            Input information about the transcripts - fields are (chromosome, strand,
            start, end, exon count, exon starts, exon ends, transcript ID, gene ID,
            gene symbol)
        include_suffix_w_tx_id : boolean
            Append parent/allele suffix to transcript names in FASTA headers of
            the output file when set to True [Default: False].

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.

        """
        #Retrieve path to the gene_files_preparation.py script.
        txptome_fasta_prep_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        txptome_fasta_prep_path = txptome_fasta_prep_path.rstrip('c')

        command = (f" python {txptome_fasta_prep_path}"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --sample_id {sample_id}"
                   f" --genome_suffix {genome_suffix}"
                   f" --genome_fasta_file_path {genome_fasta_file_path}"
                   f" --annotation_file_path {annotation_file_path}")

        if include_suffix_w_tx_id:
            command += f" --add_suffix"

        return command

    def get_validation_attributes(self, sample_id, genome_suffix, genome_fasta_file_path,
                                  annotation_file_path, include_suffix_w_tx_id=False):
        """
        Prepare attributes required by is_output_valid() function to validate
        output generated the TranscriptomeFastaPreparationStep job corresponding
        to the given input files.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to this reference genome. Used to
            construct output and log paths for this specific execution.
        genome_suffix : string
            Suffix to identify the parent/allele of the source genome. Should be
            1 or 2. This same suffix is a appended to all output files, and the
            individual transcript IDs in the transcriptome FASTA if the
            include_suffix_w_tx_id parameter is set to TRUE.
        genome_fasta_file_path : string
            Input genome fasta filename containing all the chromosomes of interest.
            No line breaks are allowed within the chromosome sequence.
        annotation_file_path : string
            Input information about the transcripts - fields are (chromosome, strand,
            start, end, exon count, exon starts, exon ends, transcript ID, gene ID,
            gene symbol)
        include_suffix_w_tx_id : boolean
            Append parent/allele suffix to transcript names in FASTA headers of
            the output file when set to True [Default: False].

        Returns
        -------
        dict
            A TranscriptomeFastaPreparationStep job's data_directory, log_directory,
            input genome_fasta_file_path, input annotation_file_path, and output
            transcriptome_fasta_file_path used when creating the transcriptome
            FASTA files.

        """
        validation_attributes = {}
        validation_attributes['data_directory'] = self.data_directory_path
        validation_attributes['log_directory'] = self.log_directory_path
        validation_attributes['sample_id'] = sample_id
        validation_attributes['genome_suffix'] = genome_suffix
        validation_attributes['genome_fasta_file_path'] = genome_fasta_file_path
        validation_attributes['annotation_file_path'] = annotation_file_path
        return validation_attributes

    @staticmethod
    def is_output_valid(validation_attributes):
        """
        Check if output of TranscriptomeFastaPreparationStep for a specific job/
        execution is correctly formed and valid, given a job's data directory,
        log directory, input genome FASTA filename, input annotation file, and
        output transcriptome FASTA filename. Prepare these attributes for a
        given jobs using the get_validation_attributes() method.

        Parameters
        ----------
        validation_attributes : dict
            A job's data_directory, log_directory, input genome_fasta_file_path,
            input annotation_file_path, and output transcriptome_fasta_file_path
            used when creating the transcriptome FASTA files.

        Returns
        -------
        boolean
            True  - TranscriptomeFastaPreparationStep output files were created
                    and are well formed.
            False - TranscriptomeFastaPreparationStep output files do not exist
                    or are missing data.

        """
        data_directory = validation_attributes['data_directory']
        log_directory_path = validation_attributes['log_directory']
        sample_id = validation_attributes['sample_id']
        genome_suffix = validation_attributes['genome_suffix']
        genome_fasta_file_path = validation_attributes['genome_fasta_file_path']
        annotation_file_path = validation_attributes['annotation_file_path']

        valid_output = False

        # Construct output filenames
        edited_genome_fasta_file_path = os.path.splitext(genome_fasta_file_path)[0] + "_edited.fa"
        trimmed_annotation_file_path = os.path.splitext(annotation_file_path)[0] + "_trimmed.txt"
        transcriptome_fasta_file_path = os.path.join(data_directory, f'sample{sample_id}',
                                                     CAMPAREE_CONSTANTS.TRANSCRIPTOME_FASTA_OUTPUT_FILENAME_PATTERN.format(genome_name=genome_suffix))
        log_file_path = os.path.join(log_directory_path, f'sample{sample_id}',
                                     CAMPAREE_CONSTANTS.TRANSCRIPTOME_FASTA_LOG_FILENAME_PATTERN.format(genome_name=genome_suffix))

        if os.path.isfile(edited_genome_fasta_file_path) and \
           os.path.isfile(trimmed_annotation_file_path) and \
           os.path.isfile(transcriptome_fasta_file_path) and \
           os.path.isfile(log_file_path):

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
        """Entry point into script when called directly.

        Parses arguments, gathers input and output filenames, and calls methods
        that perform the actual operation.

        """
        parser = argparse.ArgumentParser(description='Create transcriptome FASTA'
                                                     ' from genome sequence and annotation.')
        parser.add_argument('-l', '--log_directory_path', required=True,
                            help="Path to log directory.")
        parser.add_argument('-d', '--data_directory_path', required=True,
                            help='Path to data directory')
        parser.add_argument('--sample_id', required=True,
                            help='Sample ID associated with input genome.')
        parser.add_argument('--genome_suffix', required=True,
                            help='Suffix identifying parent/allele of source genome.')
        parser.add_argument('-g', '--genome_fasta_file_path', required=True,
                            help='Input genome FASTA filename, with no linebreaks '
                                 'within a given chromosome\'s sequence')
        parser.add_argument('-a', '--annotation_file_path', required=True,
                            help='Annotation file using coordinates matching '
                                 'the input genome FASTA file.')
        parser.add_argument('-s', '--add_suffix', action="store_true",
                            help='Append genome suffix to each transcript name '
                                 'in FASTA headers of the output file.')

        args = parser.parse_args()

        tx_fasta_prep = TranscriptomeFastaPreparationStep(args.log_directory_path,
                                                          args.data_directory_path)
        tx_fasta_prep.execute(sample_id=args.sample_id,
                              genome_suffix=args.genome_suffix,
                              genome_fasta_file_path=args.genome_fasta_file_path,
                              annotation_file_path=args.annotation_file_path,
                              include_suffix_w_tx_id=args.add_suffix)


if __name__ == "__main__":
    sys.exit(TranscriptomeFastaPreparationStep.main())
