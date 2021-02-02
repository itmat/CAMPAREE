import os
import argparse
import subprocess

from camparee.abstract_camparee_step import AbstractCampareeStep
from camparee.camparee_utils import CampareeException
from camparee.camparee_constants import CAMPAREE_CONSTANTS
from beers_utils.sample import Sample

# TODO: Add support for additional command line arguments to pass to kallisto commands.

class KallistoIndexStep(AbstractCampareeStep):
    """Wrapper around generating a kallisto transcriptome index.

    """

    KALLISTO_INDEX_DIR_PATTERN = CAMPAREE_CONSTANTS.KALLISTO_INDEX_DIR_PATTERN
    KALLISTO_INDEX_FILENAME_PATTERN = CAMPAREE_CONSTANTS.KALLISTO_INDEX_FILENAME_PATTERN
    KALLISTO_INDEX_LOG_FILENAME_PATTERN = CAMPAREE_CONSTANTS.KALLISTO_INDEX_LOG_FILENAME_PATTERN

    #The basic kallisto command used to generate transcriptome indexes.
    BASE_KALLISTO_INDEX_COMMAND = ('{kallisto_bin_path} index'
                                   ' --index={kallisto_index_file}'
                                   ' {transcriptome_fasta}')

    def __init__(self, log_directory_path, data_directory_path, parameters=None):
        """Constructor for KallistoIndexStep object.

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

    def validate(self):
        return True

    def execute(self, sample_id, genome_suffix, kallisto_bin_path, transcriptome_fasta_path):
        """Build kallisto index from the given FASTA file of transcripts.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to reference transcriptome. Used
            to construct index and log paths for this specific kallisto execution.
        genome_suffix : string
            Suffix to identify the parent/allele of the transcriptome. Should be
            1 or 2. This same suffix is a appended to all output files/directories.
        kallisto_bin_path : string
            Path to the kallisto exectuable binary.
        transcriptome_fasta_path : string
            Path to the FASTA file of transcripts, used as the basis for the
            kallisto index. This is generally prepared by the
            TranscriptomeFastaPreparationStep.

        """

        kallisto_index_dir_path = os.path.join(self.data_directory_path, f'sample{sample_id}',
                                               KallistoIndexStep.KALLISTO_INDEX_DIR_PATTERN.format(genome_name=genome_suffix))
        kallisto_index_file_path = os.path.join(self.data_directory_path, f'sample{sample_id}',
                                                KallistoIndexStep.KALLISTO_INDEX_DIR_PATTERN.format(genome_name=genome_suffix),
                                                KallistoIndexStep.KALLISTO_INDEX_FILENAME_PATTERN.format(genome_name=genome_suffix))
        log_file_path = os.path.join(self.log_directory_path, f'sample{sample_id}',
                                     KallistoIndexStep.KALLISTO_INDEX_LOG_FILENAME_PATTERN.format(genome_name=genome_suffix))

        with open(log_file_path, 'w') as log_file:

            print(f"Building kallisto indexes for transcriptome {genome_suffix} "
                           f"of sample{sample_id}.")
            log_file.write(f"Building kallisto indexes for transcriptome {genome_suffix} "
                           f"of sample{sample_id}.\n")

            log_file.write(f"Parameters:\n"
                           f"    kallisto binary path: {kallisto_bin_path}\n"
                           f"    kallisto index directory: {kallisto_index_dir_path}\n"
                           f"    kallisto index file: {kallisto_index_file_path}\n"
                           f"    input transcriptome FASTA: {transcriptome_fasta_path}\n")

            log_file.write(f"Create kallisto index directory.\n")
            if os.path.isdir(kallisto_index_dir_path):
                log_file.write(f"kallisto index directory already exists.\n")
            else:
                os.mkdir(kallisto_index_dir_path)

            kallisto_command = KallistoIndexStep.BASE_KALLISTO_INDEX_COMMAND.format(kallisto_bin_path=kallisto_bin_path,
                                                                                    kallisto_index_file=kallisto_index_file_path,
                                                                                    transcriptome_fasta=transcriptome_fasta_path)

            print(f"Running kallisto with command: {kallisto_command}")
            print(f"For full kallisto index output see {log_file_path}")
            log_file.write(f"Running kallisto with command: {kallisto_command}.\n\n")
            log_file.write("kallisto index output follows:\n")

            try:
                kallisto_result = subprocess.run(kallisto_command, shell=True, check=True,
                                                 stdout=subprocess.PIPE,
                                                 stderr=subprocess.STDOUT, #Redirect stderr to stdout.
                                                 encoding="ascii")
            except subprocess.CalledProcessError as kallisto_index_exception:
                log_file.write("\n*****ERROR: kallist index command failed:\n")
                log_file.write(f"\tExit code: {kallisto_index_exception.returncode}\n")
                log_file.write("\n*****STDOUT:\n")
                log_file.write(f"{kallisto_index_exception.stdout}\n")
                log_file.write("\n*****STDERR:\n")
                log_file.write(f"{kallisto_index_exception.stderr}\n")
                raise CampareeException(f"\nkallisto index process failed. "
                                        f"For full details see {log_file_path}\n")

            print(f"Finished generating kallisto index.\n")
            log_file.write(f"{kallisto_result.stdout}\n")
            log_file.write(f"Finished generating kallisto index.\n")
            log_file.write("ALL DONE!\n")

    def get_commandline_call(self, sample_id, genome_suffix, kallisto_bin_path, transcriptome_fasta_path):
        """
        Prepare command to execute the KallistoIndexStep from the command line,
        given all of the arugments used to run the execute() function.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to reference transcriptome. Used
            to construct index and log paths for this specific kallisto execution.
        genome_suffix : string
            Suffix to identify the parent/allele of the transcriptome. Should be
            1 or 2. This same suffix is a appended to all output files/directories.
        kallisto_bin_path : string
            Path to the kallisto binary.
        transcriptome_fasta_path : string
            Path to the FASTA file of transcripts, used as the basis for the
            kallisto index. This is generally prepared by the
            TranscriptomeFastaPreparationStep.

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.

        """

        #Retrieve path to the kallisto.py script.
        kallisto_step_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        kallisto_step_path = kallisto_step_path.rstrip('c')

        command = (f" python {kallisto_step_path} index"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --sample_id {sample_id}"
                   f" --genome_suffix {genome_suffix}"
                   f" --kallisto_bin_path {kallisto_bin_path}"
                   f" --transcriptome_fasta_file_path {transcriptome_fasta_path}")

        return command

    def get_validation_attributes(self, sample_id, genome_suffix, kallisto_bin_path, transcriptome_fasta_path):
        """
        Prepare attributes required by is_output_valid() function to validate
        output generated by the KallistoIndexStep job.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to reference transcriptome. Used
            to construct index and log paths for this specific kallisto execution.
        genome_suffix : string
            Suffix to identify the parent/allele of the transcriptome. Should be
            1 or 2. This same suffix is a appended to all output files/directories.
        kallisto_bin_path : string
            Path to the kallisto binary. [Note: this parameter is captured just
            so get_validation_attributes() accepts the same arguments as
            get_commandline_call(). It is not used here.]
        transcriptome_fasta_path : string
            Path to the FASTA file of transcripts, used as the basis for the
            kallisto index. This is generally prepared by the
            TranscriptomeFastaPreparationStep. [Note: this parameter is captured
            just so get_validation_attributes() accepts the same arguments as
            get_commandline_call(). It is not used here.]

        Returns
        -------
        dict
            A KallistoIndexStep job's data_directory, log_directory,
            corresponding sample ID, and genome_suffix.
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
        Check if output of KallistoIndexStep for a specific job/execution is
        correctly formed and valid, given a job's data directory, log directory,
        sample ID, and genome suffix. Prepare these attributes for a given job
        using the get_validation_attributes() method.

        Parameters
        ----------
        validation_attributes : dict
            A job's data_directory, log_directory, corresponding sample_id, and
            genome_suffix used when creating the kallisto index.

        Returns
        -------
        boolean
            True  - KallistoIndexStep output files were created and are well formed.
            False - KallistoIndexStep output files do not exist or are missing data.

        """

        data_directory_path = validation_attributes['data_directory']
        log_directory_path = validation_attributes['log_directory']
        sample_id = validation_attributes['sample_id']
        genome_suffix = validation_attributes['genome_suffix']

        valid_output = False

        # Construct output filenames/paths
        kallisto_index_file_path = os.path.join(data_directory_path, f'sample{sample_id}',
                                                KallistoIndexStep.KALLISTO_INDEX_DIR_PATTERN.format(genome_name=genome_suffix),
                                                KallistoIndexStep.KALLISTO_INDEX_FILENAME_PATTERN.format(genome_name=genome_suffix))
        log_file_path = os.path.join(log_directory_path, f'sample{sample_id}',
                                     KallistoIndexStep.KALLISTO_INDEX_LOG_FILENAME_PATTERN.format(genome_name=genome_suffix))

        if os.path.isfile(kallisto_index_file_path) and \
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
        """
        Entry point into class. Used when script is executed/submitted via the
        command line with the 'index' subcommand.
        """
        kallisto_index = KallistoIndexStep(log_directory_path=cmd_args.log_directory_path,
                                           data_directory_path=cmd_args.data_directory_path)
        kallisto_index.execute(sample_id=cmd_args.sample_id,
                               genome_suffix=cmd_args.genome_suffix,
                               kallisto_bin_path=cmd_args.kallisto_bin_path,
                               transcriptome_fasta_path=cmd_args.transcriptome_fasta_file_path)

class KallistoQuantStep(AbstractCampareeStep):
    """Wrapper around quantifying transript-level counts with kallisto.

    """

    KALLISTO_QUANT_DIR_PATTERN = CAMPAREE_CONSTANTS.KALLISTO_QUANT_DIR_PATTERN
    KALLISTO_ABUNDANCE_FILENAME = CAMPAREE_CONSTANTS.KALLISTO_ABUNDANCE_FILENAME
    KALLISTO_QUANT_LOG_FILENAME_PATTERN = CAMPAREE_CONSTANTS.KALLISTO_QUANT_LOG_FILENAME_PATTERN

    # The basic kallisto command used for transcript-level quants.
    BASE_KALLISTO_QUANT_COMMAND = ('{kallisto_bin_path} quant'
                                   ' --index={kallisto_index_file}'
                                   ' --output-dir={kallisto_output_dir}'
                                   ' {read_files}')

    def __init__(self, log_directory_path, data_directory_path, parameters=None):
        """Constructor for KallistoQuantStep object.

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

    def validate(self):
        return True

    def execute(self, sample, genome_suffix, kallisto_bin_path):
        """Use kallisto to generate transcript-level quantifications from fastq
        files for a given sample.

        Parameters
        ----------
        sample : Sample
            Sample containing paths for FASTQ files for quantification.
        genome_suffix : string
            Suffix to identify the parent/allele of the transcriptome. Should be
            1 or 2. This same suffix is a appended to all output files/directories.
        kallisto_bin_path : string
            Path to the kallisto exectuable binary.

        """

        kallisto_index_file_path = os.path.join(self.data_directory_path, f'sample{sample.sample_id}',
                                                KallistoIndexStep.KALLISTO_INDEX_DIR_PATTERN.format(genome_name=genome_suffix),
                                                KallistoIndexStep.KALLISTO_INDEX_FILENAME_PATTERN.format(genome_name=genome_suffix))
        kallisto_output_path = os.path.join(self.data_directory_path, f'sample{sample.sample_id}',
                                            KallistoQuantStep.KALLISTO_QUANT_DIR_PATTERN.format(genome_name=genome_suffix))
        log_file_path = os.path.join(self.log_directory_path, f'sample{sample.sample_id}',
                                     KallistoQuantStep.KALLISTO_QUANT_LOG_FILENAME_PATTERN.format(genome_name=genome_suffix))

        read_files = ' '.join(sample.fastq_file_paths)

        with open(log_file_path, 'w') as log_file:

            print(f"Running kallisto quantification for transcriptome {genome_suffix} "
                  f"of sample{sample.sample_id}")
            log_file.write(f"Running kallisto quantification for transcriptome "
                           f"{genome_suffix} of sample{sample.sample_id}.\n")

            log_file.write(f"Parameters:\n"
                           f"    kallisto binary path: {kallisto_bin_path}\n"
                           f"    kallisto index file: {kallisto_index_file_path}\n"
                           f"    kallisto output directory: {kallisto_output_path}\n"
                           f"    read files: {read_files}\n")

            log_file.write(f"Create kallisto quantification output directory.\n")
            os.mkdir(kallisto_output_path)

            kallisto_command = KallistoQuantStep.BASE_KALLISTO_QUANT_COMMAND.format(kallisto_bin_path=kallisto_bin_path,
                                                                                    kallisto_index_file=kallisto_index_file_path,
                                                                                    kallisto_output_dir=kallisto_output_path,
                                                                                    read_files=read_files)

            print(f"Running kallisto with command: {kallisto_command}")
            print(f"For full kallisto quantification output see {log_file_path}")
            log_file.write(f"Running kallisto with command: {kallisto_command}.\n\n")
            log_file.write("kallisto quantification output follows:\n")

            try:
                kallisto_result = subprocess.run(kallisto_command, shell=True, check=True,
                                                 stdout=subprocess.PIPE,
                                                 stderr=subprocess.STDOUT, #Redirect stderr to stdout.
                                                 encoding="ascii")
            except subprocess.CalledProcessError as kallisto_quant_exception:
                log_file.write("\n*****ERROR: kallist quant command failed:\n")
                log_file.write(f"\tExit code: {kallisto_quant_exception.returncode}\n")
                log_file.write("\n*****STDOUT:\n")
                log_file.write(f"{kallisto_quant_exception.stdout}\n")
                log_file.write("\n*****STDERR:\n")
                log_file.write(f"{kallisto_quant_exception.stderr}\n")
                raise CampareeException(f"\nkallisto quant process failed. "
                                        f"For full details see {log_file_path}\n")

            print(f"Finished kallisto quantification.\n")
            log_file.write(f"{kallisto_result.stdout}\n")
            log_file.write(f"Finished kallisto quantification.\n")
            log_file.write("ALL DONE!\n")

    def get_commandline_call(self, sample, genome_suffix, kallisto_bin_path):
        """
        Prepare command to execute the KallistoQuantStep from the command line,
        given all of the arugments used to run the execute() function.

        Parameters
        ----------
        sample : Sample
            Sample containing paths to FASTQ files for quantification.
        genome_suffix : string
            Suffix to identify the parent/allele of the transcriptome. Should be
            1 or 2. This same suffix is a appended to all output files/directories.
        kallisto_bin_path : string
            Path to the kallisto exectuable binary.

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.

        """

        #Retrieve path to the kallisto.py script.
        kallisto_step_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        kallisto_step_path = kallisto_step_path.rstrip('c')

        command = (f" python {kallisto_step_path} quant"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --sample '{repr(sample)}'"
                   f" --genome_suffix {genome_suffix}"
                   f" --kallisto_bin_path {kallisto_bin_path}")

        return command

    def get_validation_attributes(self, sample, genome_suffix, kallisto_bin_path):
        """
        Prepare attributes required by is_output_valid() function to validate
        output generated by the KallistoQuantStep job.

        Parameters
        ----------
        sample : Sample
            Sample containing paths to FASTQ files for quantification.
        genome_suffix : string
            Suffix to identify the parent/allele of the transcriptome. Should be
            1 or 2. This same suffix is a appended to all output files/directories.
        kallisto_bin_path : string
            Path to the kallisto exectuable binary. [Note: this parameter is
            captured just so get_validation_attributes() accepts the same
            arguments as get_commandline_call(). It is not used here.]

        Returns
        -------
        dict
            A KallistoQuantStep job's data_directory, log_directory,
            corresponding sample ID, and genome_suffix.
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
        Check if output of KallistoQuantStep for a specific job/execution is
        correctly formed and valid, given a job's data directory, log directory,
        sample ID, and genome suffix. Prepare these attributes for a given job
        using the get_validation_attributes() method.

        Parameters
        ----------
        validation_attributes : dict
            A job's data_directory, log_directory, corresponding sample_id, and
            genome_suffix used when generating transcript-level quantifications.

        Returns
        -------
        boolean
            True  - KallistoQuantStep output files were created and are well formed.
            False - KallistoQuantStep output files do not exist or are missing data.

        """

        data_directory_path = validation_attributes['data_directory']
        log_directory_path = validation_attributes['log_directory']
        sample_id = validation_attributes['sample_id']
        genome_suffix = validation_attributes['genome_suffix']

        valid_output = False

        # Construct output filenames/paths
        kallisto_output_file_path = os.path.join(data_directory_path, f'sample{sample_id}',
                                                 KallistoQuantStep.KALLISTO_QUANT_DIR_PATTERN.format(genome_name=genome_suffix),
                                                 KallistoQuantStep.KALLISTO_ABUNDANCE_FILENAME)
        log_file_path = os.path.join(log_directory_path, f'sample{sample_id}',
                                     KallistoQuantStep.KALLISTO_QUANT_LOG_FILENAME_PATTERN.format(genome_name=genome_suffix))

        if os.path.isfile(kallisto_output_file_path) and \
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
        """
        Entry point into class. Used when script is executed/submitted via the
        command line with the 'quant' subcommand.
        """
        sample = eval(cmd_args.sample)
        kallisto_quant = KallistoQuantStep(log_directory_path=cmd_args.log_directory_path,
                                           data_directory_path=cmd_args.data_directory_path)
        kallisto_quant.execute(sample=sample,
                               genome_suffix=cmd_args.genome_suffix,
                               kallisto_bin_path=cmd_args.kallisto_bin_path)

if __name__ == '__main__':
    """
    Prepare and process command line arguments. The setup belows allows for entry
    into either the KallistoIndexStep main() function or the KallistoQuantStep
    main() function based on which subcommand is specified at the command line.
    """

    parser = argparse.ArgumentParser(description='Command line wrapper around'
                                                 ' kallisto index creation and'
                                                 ' quantification.')

    subparsers = parser.add_subparsers(help="Choose one of the following:",dest="RUN_MODE", metavar="RUN_MODE")
    subparsers.required = True

    #Setup arguments for the index subcommand
    kallisto_index_subparser = subparsers.add_parser('index', help="Create kallisto index from transcriptome FASTA.",
                                                     description="Create kallisto index from transcriptome FASTA.")
    kallisto_index_subparser.set_defaults(func=KallistoIndexStep.main)
    #Send arguments for this subcommand to the KallistoIndexStep's main() method.
    required_named_kallisto_index_subparser = kallisto_index_subparser.add_argument_group('Required named arguments')
    required_named_kallisto_index_subparser.add_argument('-l', '--log_directory_path', required=True,
                                                         help='Directory in which to save logging files.')
    required_named_kallisto_index_subparser.add_argument('-d', '--data_directory_path', required=True,
                                                         help='Directory in which to save output files.')
    required_named_kallisto_index_subparser.add_argument('--sample_id', required=True,
                                                         help='Sample ID associated with input genome.')
    required_named_kallisto_index_subparser.add_argument('--genome_suffix', required=True,
                                                         help='Suffix identifying parent/allele of source genome.')
    required_named_kallisto_index_subparser.add_argument('--kallisto_bin_path', required=True,
                                                         help='Full path to kallisto executable binary')
    required_named_kallisto_index_subparser.add_argument('--transcriptome_fasta_file_path', required=True,
                                                         help='Input transcriptome in FASTA format.')

    #Setup arguments from the quantification subcommand
    kallisto_quant_subparser = subparsers.add_parser('quant', help="Run kallisto transcript-level quantification.",
                                                   description="Run kallisto transcript-level quantification.")
    #Send arguments for this subcommand to the KallistoQuantStep's main() method.
    kallisto_quant_subparser.set_defaults(func=KallistoQuantStep.main)
    required_named_kallisto_quant_subparser = kallisto_quant_subparser.add_argument_group('Required named arguments')
    required_named_kallisto_quant_subparser.add_argument('-l', '--log_directory_path', required=True,
                                                         help='Directory in which to save logging files.')
    required_named_kallisto_quant_subparser.add_argument('-d', '--data_directory_path', required=True,
                                                         help='Directory in which to save output files.')
    required_named_kallisto_quant_subparser.add_argument('--sample', required=True,
                                                         help='String representation of a Sample object.')
    required_named_kallisto_quant_subparser.add_argument('--genome_suffix', required=True,
                                                         help='Suffix identifying parent/allele of source genome.')
    required_named_kallisto_quant_subparser.add_argument('--kallisto_bin_path', required=True,
                                                         help='Full path to kallisto executable binary')

    args = parser.parse_args()
    args.func(args)
