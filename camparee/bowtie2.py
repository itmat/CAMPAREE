import os
import sys
import argparse
import subprocess
import json

from camparee.abstract_camparee_step import AbstractCampareeStep
from camparee.camparee_utils import CampareeException
from camparee.camparee_constants import CAMPAREE_CONSTANTS
from beers_utils.sample import Sample

# TODO: Add support for additional command line arguments to pass to Bowtie2 commands.

class Bowtie2IndexStep(AbstractCampareeStep):
    """Wrapper around generating a Bowtie2 index.

    """

    BOWTIE2_INDEX_DIR_PATTERN = CAMPAREE_CONSTANTS.BOWTIE2_INDEX_DIR_PATTERN
    BOWTIE2_INDEX_PREFIX_PATTERN = CAMPAREE_CONSTANTS.BOWTIE2_INDEX_PREFIX_PATTERN
    BOWTIE2_INDEX_LOG_FILENAME_PATTERN = CAMPAREE_CONSTANTS.BOWTIE2_INDEX_LOG_FILENAME_PATTERN

    #The basic Bowtie2 command used to generate indexes from a given FASTA.
    BASE_BOWTIE2_INDEX_COMMAND = ('{bowtie2_bin_dir}/bowtie2-build'
                                  ' --threads {num_bowtie2_threads}'
                                  ' {bowtie2_cmd_options}'
                                  ' {reference_fasta}'
                                  ' {output_index_prefix}')

    def __init__(self, log_directory_path, data_directory_path, parameters=dict()):
        """Constructor for Bowtie2IndexStep object.

        Parameters
        ----------
        data_directory_path: string
            Full path to data directory
        log_directory_path : string
            Full path to log directory.
        parameters : dict
            [Optional] Dictionary of Bowtie2 parameters specified by the config
            file (Note, the "num_bowtie_threads" entry in the config file maps
            to the bowtie2 "--threads" command line parameter).

        """
        self.data_directory_path = data_directory_path
        self.log_directory_path = log_directory_path
        self.num_bowtie2_threads = parameters.pop('num_bowtie_threads', 1)
        # Remaining parameters (if any) aside from  "num_bowtie_threads"
        self.bowtie2_cmd_options = parameters

    def validate(self):
        # The value given to the "--threads" parameter is specified by the
        # num_bowtie_threads entry in the config file.
        invalid_bowtie2_parameters = ["--threads"]
        for key, value in self.bowtie2_cmd_options.items():
            if not key.startswith("-"):
                print(f"Bowtie2 index parameter {key} with value {value} needs"
                      f" to be a Bowtie2 option starting with single (-) or double"
                      f" dashes (--).",
                      sys.stderr)
                return False
            if key in invalid_bowtie2_parameters:
                print(f"Bowtie2 index parameter {key} with value {value} cannot"
                      f" be used as a Bowtie2 option since the value is either"
                      f" hard-coded by this script, or explicitly specfied"
                      f" elsewhere in the config file.")
                return False

        return True

    def execute(self, sample_id, genome_suffix, bowtie2_bin_dir, transcriptome_fasta_path):
        """Build Bowtie2 index from the given FASTA file of transcripts.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to reference transcriptome. Used
            to construct index and log paths for this specific Bowtie2 execution.
        genome_suffix : string
            Suffix to identify the parent/allele of the transcriptome. Should be
            1 or 2. This same suffix is a appended to all output files/directories.
        bowtie2_bin_dir : string
            Path to the directory containing the bowtie2-build exectuable.
        transcriptome_fasta_path : string
            Path to the FASTA file of transcripts, used as the basis for the
            Bowtie2 index. This is generally prepared by the
            TranscriptomeFastaPreparationStep.

        """
        bowtie2_index_dir_path = os.path.join(self.data_directory_path, f'sample{sample_id}',
                                               Bowtie2IndexStep.BOWTIE2_INDEX_DIR_PATTERN.format(genome_name=genome_suffix))
        bowtie2_index_file_prefix = os.path.join(self.data_directory_path, f'sample{sample_id}',
                                                 Bowtie2IndexStep.BOWTIE2_INDEX_DIR_PATTERN.format(genome_name=genome_suffix),
                                                 Bowtie2IndexStep.BOWTIE2_INDEX_PREFIX_PATTERN.format(genome_name=genome_suffix))
        log_file_path = os.path.join(self.log_directory_path, f'sample{sample_id}',
                                     Bowtie2IndexStep.BOWTIE2_INDEX_LOG_FILENAME_PATTERN.format(genome_name=genome_suffix))

        with open(log_file_path, 'w') as log_file:

            print(f"Building Bowtie2 indexes for transcriptome {genome_suffix} "
                           f"of sample{sample_id}.")
            log_file.write(f"Building Bowtie2 indexes for transcriptome {genome_suffix} "
                           f"of sample{sample_id}.\n")

            log_file.write(f"Parameters:\n"
                           f"    Bowtie2 binary directory: {bowtie2_bin_dir}\n"
                           f"    Bowtie2 index directory: {bowtie2_index_dir_path}\n"
                           f"    Bowtie2 index file prefix: {bowtie2_index_file_prefix}\n"
                           f"    Input transcriptome FASTA: {transcriptome_fasta_path}\n"
                           f"    Number of Bowtie2 threads: {self.num_bowtie2_threads}\n")

            log_file.write(f"Create Bowtie2 index directory.\n")
            os.mkdir(bowtie2_index_dir_path)

            bwt2_cmd_options = ' '.join( f"{key} {value}" for key,value in self.bowtie2_cmd_options.items() )

            bowtie2_command = Bowtie2IndexStep.BASE_BOWTIE2_INDEX_COMMAND.format(bowtie2_bin_dir=bowtie2_bin_dir,
                                                                                 num_bowtie2_threads=self.num_bowtie2_threads,
                                                                                 bowtie2_cmd_options=bwt2_cmd_options,
                                                                                 output_index_prefix=bowtie2_index_file_prefix,
                                                                                 reference_fasta=transcriptome_fasta_path)

            print(f"Running Bowtie2 with command: {bowtie2_command}")
            print("Bowtie2 index output follows:\n")
            log_file.write(f"Running Bowtie2 with command: {bowtie2_command}.\n\n")
            log_file.write("Bowtie2 index output follows:\n")

            try:
                bowtie2_result = subprocess.run(bowtie2_command, shell=True, check=True,
                                                stdout=subprocess.PIPE,
                                                stderr=subprocess.STDOUT, #Redirect stderr to stdout.
                                                encoding="ascii")
            except subprocess.CalledProcessError as bowtie2_index_exception:
                log_file.write("\n*****ERROR: Bowtie2 index command failed:\n")
                log_file.write(f"\tExit code: {bowtie2_index_exception.returncode}\n")
                log_file.write("\n*****STDOUT:\n")
                log_file.write(f"{bowtie2_index_exception.stdout}\n")
                log_file.write("\n*****STDERR:\n")
                log_file.write(f"{bowtie2_index_exception.stderr}\n")
                raise CampareeException(f"\nBowtie2 index process failed. "
                                        f"For full details see {log_file_path}\n")

            print(bowtie2_result.stdout)
            print(f"\nFinished generating Bowtie2 index.\n")
            log_file.write(f"{bowtie2_result.stdout}\n")
            log_file.write(f"\nFinished generating Bowtie2 index.\n")
            log_file.write("ALL DONE!\n")

    def get_commandline_call(self, sample_id, genome_suffix, bowtie2_bin_dir, transcriptome_fasta_path):
        """Prepare command to execute the Bowtie2IndexStep from the command line,
        given all of the arugments used to run the execute() function.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to reference transcriptome. Used
            to construct index and log paths for this specific Bowtie2 execution.
        genome_suffix : string
            Suffix to identify the parent/allele of the transcriptome. Should be
            1 or 2. This same suffix is a appended to all output files/directories.
        bowtie2_bin_dir : string
            Path to the directory containing the bowtie2-build exectuable.
        transcriptome_fasta_path : string
            Path to the FASTA file of transcripts, used as the basis for the
            Bowtie2 index. This is generally prepared by the
            TranscriptomeFastaPreparationStep.

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.

        """

        #Retrieve path to the bowtie2.py script.
        bowtie2_step_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        bowtie2_step_path = bowtie2_step_path.rstrip('c')

        # TODO: Explore alternative to json for passing dictionary via command line.
        #       Eval could be dangerous for this, since the user has complete control
        #       over what gets entered as a bowtie2 parameter through the config file.

        command = (f" python {bowtie2_step_path} index"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --sample_id {sample_id}"
                   f" --genome_suffix {genome_suffix}"
                   f" --bowtie2_bin_dir {bowtie2_bin_dir}"
                   f" --transcriptome_fasta_file_path {transcriptome_fasta_path}"
                   f" --num_bowtie2_threads {self.num_bowtie2_threads}"
                   f" --bowtie2_parameters '{json.dumps(self.bowtie2_cmd_options)}'")
        return command

    def get_validation_attributes(self, sample_id, genome_suffix, bowtie2_bin_dir, transcriptome_fasta_path):
        """
        Prepare attributes required by is_output_valid() function to validate
        output generated by the Bowtie2IndexStep job.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to reference transcriptome. Used
            to construct index and log paths for this specific Bowtie2 execution.
        genome_suffix : string
            Suffix to identify the parent/allele of the transcriptome. Should be
            1 or 2. This same suffix is a appended to all output files/directories.
        bowtie2_bin_dir : string
            Path to the directory containing the bowtie2-build exectuable. [Note:
            this parameter is captured just so get_validation_attributes() accepts
            the same arguments as get_commandline_call(). It is not used here.]
        transcriptome_fasta_path : string
            Path to the FASTA file of transcripts, used as the basis for the
            Bowtie2 index. This is generally prepared by the
            TranscriptomeFastaPreparationStep. [Note: this parameter is captured
            just so get_validation_attributes() accepts the same arguments as
            get_commandline_call(). It is not used here.]

        Returns
        -------
        dict
            A Bowtie2IndexStep job's data_directory, log_directory, corresponding
            sample ID, and genome_suffix.
        """
        validation_attributes = {}
        validation_attributes['data_directory'] = self.data_directory_path
        validation_attributes['log_directory'] = self.log_directory_path
        validation_attributes['sample_id'] = sample_id
        validation_attributes['genome_suffix'] = genome_suffix
        return validation_attributes

    @staticmethod
    def is_output_valid(validation_attributes):
        """
        Check if output of Bowtie2IndexStep for a specific job/execution is
        correctly formed and valid, given a job's data directory, log directory,
        sample ID, and genome suffix. Prepare these attributes for a given job
        using the get_validation_attributes() method.

        Parameters
        ----------
        validation_attributes : dict
            A job's data_directory, log_directory, corresponding sample_id, and
            genome_suffix used when creating the Bowtie2 index.

        Returns
        -------
        boolean
            True  - Bowtie2IndexStep output files were created and are well formed.
            False - Bowtie2IndexStep output files do not exist or are missing data.

        """

        data_directory_path = validation_attributes['data_directory']
        log_directory_path = validation_attributes['log_directory']
        sample_id = validation_attributes['sample_id']
        genome_suffix = validation_attributes['genome_suffix']

        valid_output = False

        # Construct output filenames/paths
        log_file_path = os.path.join(log_directory_path, f'sample{sample_id}',
                                     Bowtie2IndexStep.BOWTIE2_INDEX_LOG_FILENAME_PATTERN.format(genome_name=genome_suffix))
        bowtie2_index_file_prefix = os.path.join(data_directory_path, f'sample{sample_id}',
                                                 Bowtie2IndexStep.BOWTIE2_INDEX_DIR_PATTERN.format(genome_name=genome_suffix),
                                                 Bowtie2IndexStep.BOWTIE2_INDEX_PREFIX_PATTERN.format(genome_name=genome_suffix))

        # TODO: Identify index files are missing in the event of a failed validation.

        # Note, bowtie2-build should produce 6 different index files. They all
        # should exist.
        if os.path.isfile(bowtie2_index_file_prefix + ".1.bt2") and \
           os.path.isfile(bowtie2_index_file_prefix + ".2.bt2") and \
           os.path.isfile(bowtie2_index_file_prefix + ".3.bt2") and \
           os.path.isfile(bowtie2_index_file_prefix + ".4.bt2") and \
           os.path.isfile(bowtie2_index_file_prefix + ".rev.1.bt2") and \
           os.path.isfile(bowtie2_index_file_prefix + ".rev.2.bt2") and \
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
    def main(cmd_args):
        """Entry point into class. Used when script is executed/submitted via
        the command line with the 'index' subcommand.
        """
        parameters = json.loads(cmd_args.bowtie2_parameters)
        bowtie2_index = Bowtie2IndexStep(log_directory_path=cmd_args.log_directory_path,
                                         data_directory_path=cmd_args.data_directory_path,
                                         parameters=parameters)
        bowtie2_index.execute(sample_id=cmd_args.sample_id,
                              genome_suffix=cmd_args.genome_suffix,
                              bowtie2_bin_dir=cmd_args.bowtie2_bin_dir,
                              transcriptome_fasta_path=cmd_args.transcriptome_fasta_file_path)

class Bowtie2AlignStep(AbstractCampareeStep):
    """Wrapper around aligning reads with Bowtie2

    """

    BOWTIE2_ALIGN_FILENAME_PATTERN = CAMPAREE_CONSTANTS.BOWTIE2_ALIGN_FILENAME_PATTERN
    BOWTIE2_ALIGN_LOG_FILENAME_PATTERN = CAMPAREE_CONSTANTS.BOWTIE2_ALIGN_LOG_FILENAME_PATTERN

    # TODO: Update this script to gracefully handle both one and two FASTQ files
    #       for input (currently only works with two FASTQ files).

    #The basic Bowtie2 command used to generate indexes from a given FASTA.
    BASE_BOWTIE2_ALIGN_COMMAND = ('{bowtie2_bin_dir}/bowtie2'
                                  ' --very-sensitive'
                                  ' --threads {num_bowtie2_threads}'
                                  ' {bowtie2_cmd_options}'
                                  ' -x {bowtie2_index_prefix}'
                                  ' -1 {first_read_fastq}'
                                  ' -2 {second_read_fastq}'
                                  ' -S {output_sam_file}')

    def __init__(self, log_directory_path, data_directory_path, parameters=dict()):
        """Constructor for Bowtie2AlignStep object.

        Parameters
        ----------
        data_directory_path: string
           Full path to data directory
        log_directory_path : string
           Full path to log directory.
        parameters : dict
            [Optional] Dictionary of Bowtie2 parameters specified by the config
            file (Note, the "num_bowtie_threads" entry in the config file maps
            to the bowtie2 "--threads" command line parameter).

        """
        self.data_directory_path = data_directory_path
        self.log_directory_path = log_directory_path
        self.num_bowtie2_threads = parameters.pop('num_bowtie_threads', 1)
        # Remaining parameters (if any) aside from  "num_bowtie_threads"
        self.bowtie2_cmd_options = parameters

    def validate(self):
        """Check all given Bowtie2 parameters are correctly formed (i.e. start
        with single or double dash), and do not conflict with any that are
        explicitly specified by this script (--very-sensitive, -x, -1, -2, -S),
        or elsewhere in the config file (--threads).

        """
        # These are parameters this script specifies directly. Most of these are
        # for specifying the index, input fastq(s), and output SAM filename.
        invalid_bowtie2_parameters = ["--very-sensitive", "-x", "-1", "-2", "-S", "--threads"]
        for key, value in self.bowtie2_cmd_options.items():
            if not key.startswith("-"):
                print(f"Bowtie2 align parameter {key} with value {value} needs"
                      f" to be a Bowtie2 option starting with single (-) or double"
                      f" dashes (--).",
                      sys.stderr)
                return False
            if key in invalid_bowtie2_parameters:
                print(f"Bowtie2 index parameter {key} with value {value} cannot"
                      f" be used as a Bowtie2 option since the value is either"
                      f" hard-coded by this script, or explicitly specfied"
                      f" elsewhere in the config file.")
                return False

        return True

    def execute(self, sample, genome_suffix, bowtie2_bin_dir):
        """Use Bowtie2 to align fastq files for a given sample to the refrence
        transcriptome.

        Parameters
        ----------
        sample : Sample
            Sample containing paths for FASTQ files for alignment.
        genome_suffix : string
            Suffix to identify the parent/allele of the transcriptome. Should be
            1 or 2. This same suffix is a appended to all output files/directories.
        bowtie2_bin_dir : string
            Path to the directory containing the bowtie2 exectuable.

        """

        bowtie2_index_file_prefix = os.path.join(self.data_directory_path, f'sample{sample.sample_id}',
                                                 Bowtie2IndexStep.BOWTIE2_INDEX_DIR_PATTERN.format(genome_name=genome_suffix),
                                                 Bowtie2IndexStep.BOWTIE2_INDEX_PREFIX_PATTERN.format(genome_name=genome_suffix))
        bowtie2_output_file_path = os.path.join(self.data_directory_path, f'sample{sample.sample_id}',
                                                Bowtie2AlignStep.BOWTIE2_ALIGN_FILENAME_PATTERN.format(genome_name=genome_suffix))
        log_file_path = os.path.join(self.log_directory_path, f'sample{sample.sample_id}',
                                     Bowtie2AlignStep.BOWTIE2_ALIGN_LOG_FILENAME_PATTERN.format(genome_name=genome_suffix))

        fastq_file_1, fastq_file_2 = sample.fastq_file_paths

        with open(log_file_path, 'w') as log_file:

            print(f"Running Bowtie2 alignment to transcriptome {genome_suffix} "
                  f"of sample{sample.sample_id}")
            log_file.write(f"Running Bowtie2 alignment to transcriptome "
                           f"{genome_suffix} of sample{sample.sample_id}.\n")

            log_file.write(f"Parameters:\n"
                           f"    Bowtie2 binary directory: {bowtie2_bin_dir}\n"
                           f"    Bowtie2 index file prefix: {bowtie2_index_file_prefix}\n"
                           f"    Bowtie2 output SAM file: {bowtie2_output_file_path}\n"
                           f"    Read 1 FASTQ: {fastq_file_1}\n"
                           f"    Read 2 FASTQ: {fastq_file_2}\n"
                           f"    Number of Bowtie2 threads: {self.num_bowtie2_threads}\n")

            bwt2_cmd_options = ' '.join( f"{key} {value}" for key,value in self.bowtie2_cmd_options.items() )

            bowtie2_command = Bowtie2AlignStep.BASE_BOWTIE2_ALIGN_COMMAND.format(bowtie2_bin_dir=bowtie2_bin_dir,
                                                                                 num_bowtie2_threads=self.num_bowtie2_threads,
                                                                                 bowtie2_cmd_options=bwt2_cmd_options,
                                                                                 bowtie2_index_prefix=bowtie2_index_file_prefix,
                                                                                 first_read_fastq=fastq_file_1,
                                                                                 second_read_fastq=fastq_file_2,
                                                                                 output_sam_file=bowtie2_output_file_path)

            print(f"Running Bowtie2 with command: {bowtie2_command}")
            print("Bowtie2 alignment output follows:\n")
            log_file.write(f"Running Bowtie2 with command: {bowtie2_command}.\n\n")
            log_file.write("Bowtie2 alignment output follows:\n")

            try:
                bowtie2_result = subprocess.run(bowtie2_command, shell=True, check=True,
                                                stdout=subprocess.PIPE,
                                                stderr=subprocess.STDOUT, #Redirect stderr to stdout.
                                                encoding="ascii")
            except subprocess.CalledProcessError as bowtie2_align_exception:
                log_file.write("\n*****ERROR: Bowtie2 alignment command failed:\n")
                log_file.write(f"\tExit code: {bowtie2_align_exception.returncode}\n")
                log_file.write("\n*****STDOUT:\n")
                log_file.write(f"{bowtie2_align_exception.stdout}\n")
                log_file.write("\n*****STDERR:\n")
                log_file.write(f"{bowtie2_align_exception.stderr}\n")
                raise CampareeException(f"\nBowtie2 alignment process failed. "
                                        f"For full details see {log_file_path}\n")

            print(bowtie2_result.stdout)
            print(f"\nFinished Bowtie2 alignment.\n")
            log_file.write(f"{bowtie2_result.stdout}\n")
            log_file.write(f"\nFinished Bowtie2 alignment.\n")
            log_file.write("ALL DONE!\n")

    def get_commandline_call(self, sample, genome_suffix, bowtie2_bin_dir):
        """Prepare command to execute the Bowtie2AlignStep from the command line,
        given all of the arugments used to run the execute() function.

        Parameters
        ----------
        sample : Sample
            Sample containing paths for FASTQ files for alignment.
        genome_suffix : string
            Suffix to identify the parent/allele of the transcriptome. Should be
            1 or 2. This same suffix is a appended to all output files/directories.
        bowtie2_bin_dir : string
            Path to the directory containing the bowtie2 exectuable.

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.

        """

        #Retrieve path to the bowtie2.py script.
        bowtie2_step_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        bowtie2_step_path = bowtie2_step_path.rstrip('c')

        # TODO: Explore alternative to json for passing dictionary via command line.
        #       Eval could be dangerous for this, since the user has complete control
        #       over what gets entered as a bowtie2 parameter through the config file.

        command = (f" python {bowtie2_step_path} align"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --sample '{repr(sample)}'"
                   f" --genome_suffix {genome_suffix}"
                   f" --bowtie2_bin_dir {bowtie2_bin_dir}"
                   f" --num_bowtie2_threads {self.num_bowtie2_threads}"
                   f" --bowtie2_parameters '{json.dumps(self.bowtie2_cmd_options)}'")

        return command

    def get_validation_attributes(self, sample, genome_suffix, bowtie2_bin_dir):
        """
        Prepare attributes required by is_output_valid() function to validate
        output generated by the Bowtie2AlignStep job.

        Parameters
        ----------
        sample : Sample
            Sample containing paths for FASTQ files for alignment. [Note: only
            the sample_id is used, but the full Sample object is required here
            so get_validation_attributes() accepts the same arguments as
            get_commandline_call().]
        genome_suffix : string
            Suffix to identify the parent/allele of the transcriptome. Should be
            1 or 2. This same suffix is a appended to all output files/directories.
        bowtie2_bin_dir : string
            Path to the directory containing the bowtie2 exectuable. [Note: this
            parameter is captured just so get_validation_attributes() accepts
            the same arguments as get_commandline_call(). It is not used here.]

        Returns
        -------
        dict
            A Bowtie2AlignStep job's data_directory, log_directory, corresponding
            sample ID, and genome_suffix.
        """
        validation_attributes = {}
        validation_attributes['data_directory'] = self.data_directory_path
        validation_attributes['log_directory'] = self.log_directory_path
        validation_attributes['sample_id'] = sample.sample_id
        validation_attributes['genome_suffix'] = genome_suffix
        return validation_attributes

    @staticmethod
    def is_output_valid(validation_attributes):
        """
        Check if output of Bowtie2AlignStep for a specific job/execution is
        correctly formed and valid, given a job's data directory, log directory,
        sample ID, and genome suffix. Prepare these attributes for a given job
        using the get_validation_attributes() method.

        Parameters
        ----------
        validation_attributes : dict
            A job's data_directory, log_directory, corresponding sample_id, and
            genome_suffix used when aligning reads with Bowtie2.

        Returns
        -------
        boolean
            True  - Bowtie2AlignStep output files were created and are well formed.
            False - Bowtie2AlignStep output files do not exist or are missing data.

        """

        data_directory_path = validation_attributes['data_directory']
        log_directory_path = validation_attributes['log_directory']
        sample_id = validation_attributes['sample_id']
        genome_suffix = validation_attributes['genome_suffix']

        valid_output = False

        # Construct output filenames/paths
        bowtie2_output_file_path = os.path.join(data_directory_path, f'sample{sample_id}',
                                                Bowtie2AlignStep.BOWTIE2_ALIGN_FILENAME_PATTERN.format(genome_name=genome_suffix))
        log_file_path = os.path.join(log_directory_path, f'sample{sample_id}',
                                     Bowtie2AlignStep.BOWTIE2_ALIGN_LOG_FILENAME_PATTERN.format(genome_name=genome_suffix))

        if os.path.isfile(bowtie2_output_file_path) and \
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
    def main(cmd_args):
        """Entry point into class. Used when script is executed/submitted via
        the command line with the 'align' subcommand.
        """
        sample = eval(cmd_args.sample)
        parameters = json.loads(cmd_args.bowtie2_parameters)
        parameters['num_bowtie_threads'] = cmd_args.num_bowtie2_threads
        bowtie2_align = Bowtie2AlignStep(log_directory_path=cmd_args.log_directory_path,
                                         data_directory_path=cmd_args.data_directory_path,
                                         parameters=parameters)
        bowtie2_align.execute(sample=sample,
                              genome_suffix=cmd_args.genome_suffix,
                              bowtie2_bin_dir=cmd_args.bowtie2_bin_dir)

if __name__ == '__main__':
    """
    Prepare and process command line arguments. The setup below allows for entry
    into either the Bowtie2IndexStep main() method or the Bowtie2AlignStep
    main() method based on which subcommand is specified at the command line.
    """

    parser = argparse.ArgumentParser(description='Command line wrapper around'
                                                 ' Bowtie2 index creation and'
                                                 ' alignment.')

    subparsers = parser.add_subparsers(help="Choose one of the following:",dest="RUN_MODE", metavar="RUN_MODE")
    subparsers.required = True

    #Setup arguments for the index subcommand
    bowtie2_index_subparser = subparsers.add_parser('index', help="Create Bowtie2 index from transcriptome FASTA.",
                                                     description="Create Bowtie2 index from transcriptome FASTA.")
    bowtie2_index_subparser.set_defaults(func=Bowtie2IndexStep.main)
    #Send arguments for this subcommand to the Bowtie2IndexStep's main() method.
    required_named_bowtie2_index_subparser = bowtie2_index_subparser.add_argument_group('Required named arguments')
    required_named_bowtie2_index_subparser.add_argument('-l', '--log_directory_path', required=True,
                                                        help='Directory in which to save logging files.')
    required_named_bowtie2_index_subparser.add_argument('-d', '--data_directory_path', required=True,
                                                        help='Directory in which to save output files.')
    required_named_bowtie2_index_subparser.add_argument('--sample_id', required=True,
                                                        help='Sample ID associated with input genome.')
    required_named_bowtie2_index_subparser.add_argument('--genome_suffix', required=True,
                                                        help='Suffix identifying parent/allele of source genome.')
    required_named_bowtie2_index_subparser.add_argument('--bowtie2_bin_dir', required=True,
                                                        help='Full path to directory containing bowtie2-build '
                                                             'executable.')
    required_named_bowtie2_index_subparser.add_argument('--transcriptome_fasta_file_path', required=True,
                                                        help='Input transcriptome in FASTA format.')
    required_named_bowtie2_index_subparser.add_argument('--num_bowtie2_threads', type=int, default=1, required=False,
                                                        help='Number of threads to use when running Bowtie2,')
    required_named_bowtie2_index_subparser.add_argument('--bowtie2_parameters', required=False,
                                                        help="Jsonified Bowtie2 index parameters (excluding "
                                                             "--threads).")

    #Setup arguments from the alignment subcommand
    bowtie2_align_subparser = subparsers.add_parser('align', help="Run Bowtie2 alignment to transcriptome.",
                                                   description="Run Bowtie2 alignment to transcriptome.")
    #Send arguments for this subcommand to the Bowtie2AlignStep's main() method.
    bowtie2_align_subparser.set_defaults(func=Bowtie2AlignStep.main)
    required_named_bowtie2_align_subparser = bowtie2_align_subparser.add_argument_group('Required named arguments')
    required_named_bowtie2_align_subparser.add_argument('-l', '--log_directory_path', required=True,
                                                        help='Directory in which to save logging files.')
    required_named_bowtie2_align_subparser.add_argument('-d', '--data_directory_path', required=True,
                                                        help='Directory in which to save output files.')
    required_named_bowtie2_align_subparser.add_argument('--sample', required=True,
                                                        help='String representation of a Sample object.')
    required_named_bowtie2_align_subparser.add_argument('--genome_suffix', required=True,
                                                        help='Suffix identifying parent/allele of source genome.')
    required_named_bowtie2_align_subparser.add_argument('--bowtie2_bin_dir', required=True,
                                                        help='Full path to directory containing bowtie2 '
                                                             'executable.')
    required_named_bowtie2_align_subparser.add_argument('--num_bowtie2_threads', type=int, default=1, required=False,
                                                        help='Number of threads to use when running Bowtie2,')
    required_named_bowtie2_align_subparser.add_argument('--bowtie2_parameters', required=False,
                                                        help="Jsonified Bowtie2 index parameters (excluding "
                                                             "--threads).")

    args = parser.parse_args()
    args.func(args)
