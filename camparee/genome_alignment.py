import os
import sys
import argparse
import shutil
import subprocess
import re #Probably won't need this if we switch to an LSF API
import json

import pysam

from beers_utils.sample import Sample
from camparee.abstract_camparee_step import AbstractCampareeStep
from camparee.camparee_utils import CampareeException
from camparee.camparee_constants import CAMPAREE_CONSTANTS

class GenomeAlignmentStep(AbstractCampareeStep):

    #The basic STAR command used to align the given fastq files the chosen genome.
    BASE_STAR_COMMAND = ('{star_bin_path}'
                         ' --outFileNamePrefix {out_file_prefix}'
                         ' --genomeDir {star_index_path}'
                         ' --runMode alignReads'
                         ' --outSAMtype BAM SortedByCoordinate'
                         ' {star_cmd_options}'
                         ' --readFilesIn {read_files}')

    def __init__(self, log_directory_path, data_directory_path, parameters = dict()):
        self.log_directory_path = log_directory_path
        self.data_directory_path = data_directory_path
        self.star_cmd_options = parameters

    def validate(self):
        invalid_STAR_parameters = ["--outFileNamePrefix", "--genomeDir", "--runMode", "--outSAMtype", "--readFilesIn"]
        for key, value in self.star_cmd_options.items():
            if not key.startswith("--"):
                print(f"Genome Alignment parameter {key} with value {value} needs to be"
                      f"a STAR option starting with double dashes --", sys.stderr)
                return False
            if key in invalid_STAR_parameters:
                print(f"Genome Alignment parameter {key} with value {value} cannot be"
                      f" used as a STAR option since the value is already determined by the genome alignment script.")
                return False

        return True

    def execute(self, sample, star_index_directory_path, star_bin_path):
        """
        Use STAR to align fastq files for a given sample, to the reference genome.

        Parameters
        ----------
        sample : Sample
            Sample containing paths for FASTQ files to align, or pre-aligned
            BAM file.
        star_index_directory_path : string
            Path to directory containing STAR index.
        star_bin_path : string
            Path to STAR executable binary.

        """

        # Check if user provided bam files, so we don't re-run the alignment
        if sample.bam_file_path:
            if os.path.isfile(sample.bam_file_path):
                print(f"sample{sample.sample_id} {sample.sample_name} is already aligned.")
                return
            else:
                raise CampareeException(f"BAM file provided for sample{sample.sample_id} - {sample.sample_name}"
                                        f" does note exist at the following path:\n"
                                        f"{sample.bam_file_path}")


        read_files = ' '.join(sample.fastq_file_paths)
        out_file_prefix = os.path.join(self.data_directory_path, f"sample{sample.sample_id}",
                                       CAMPAREE_CONSTANTS.DEFAULT_STAR_OUTPUT_PREFIX)
        star_cmd_options = ' '.join( f"{key} {value}" for key,value in self.star_cmd_options.items() )

        star_command = self.BASE_STAR_COMMAND.format(star_bin_path=star_bin_path,
                                                     out_file_prefix=out_file_prefix,
                                                     star_index_path=star_index_directory_path,
                                                     star_cmd_options=star_cmd_options,
                                                     read_files=read_files)

        print(f"Starting STAR on sample {sample.sample_name}.")
        star_result = subprocess.run(star_command, shell=True, check=True)
        print(f"Finished running STAR on {sample.sample_name}.")

    def get_genome_bam_path(self, sample):
        """
        Determine whether user provided a BAM file for the given sample, and
        return either this path, or the default path used by the GenomeAlignment
        step.

        Parameters
        ----------
        sample : Sample
            Sample containing paths for FASTQ files to align, or pre-aligned
            BAM file.

        Returns
        -------
        string
            Path to BAM file associated with this sample. Either path given by
            user or default path used by GenomeAlignment step.

        """

        if sample.bam_file_path:
            if os.path.isfile(sample.bam_file_path):
                return sample.bam_file_path
            else:
                raise CampareeException(f"BAM file provided for sample{sample.sample_id} - {sample.sample_name}"
                                        f" does note exist at the following path:\n"
                                        f"{sample.bam_file_path}")
        else:
            default_bam_path = os.path.join(self.data_directory_path,
                                            f"sample{sample.sample_id}",
                                            f"{CAMPAREE_CONSTANTS.DEFAULT_STAR_BAM_FILENAME}")
            return default_bam_path

    def get_commandline_call(self, sample, star_index_directory_path, star_bin_path):
        """
        Prepare command to execute the GenomeAlignment from the command line, given
        all of the arugments used to run the execute() function.

        Parameters
        ----------
        sample : Sample
            Sample containing paths for FASTQ files to align, or pre-aligned
            BAM file.
        star_index_directory_path : string
            Path to directory containing STAR index.
        star_bin_path : string
            Path to STAR executable binary.

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.

        """

        #Retrieve path to the genome_alignment.py script.
        genome_align_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        genome_align_path = genome_align_path.rstrip('c')

        command = (f" python {genome_align_path} align"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --star_index_directory_path {star_index_directory_path}"
                   f" --sample '{repr(sample)}'"
                   f" --star_bin_path {star_bin_path}"
                   f" --star_parameters '{json.dumps(self.star_cmd_options)}'")

        return command

    def get_validation_attributes(self, sample, star_index_directory_path, star_bin_path):
        """
        Prepare attributes required by is_output_valid() function to validate
        output generated the STAR genome job corresponding to the given sample.

        Parameters
        ----------
        sample : Sample
            Sample defining the FASTQ files to be aligned, or the pre-aligned BAM.
        star_index_directory_path : string
            Path to directory containing STAR index. [Note: this parameter is
            captured just so get_validation_attributes() accepts the same arguments
            as get_commandline_call(). It is not used here.]
        star_bin_path : string
            Path to STAR executable binary. [Note: this parameter is captured
            just so get_validation_attributes() accepts the same arguments as
            get_commandline_call(). It is not used here.]

        Returns
        -------
        dict
            The GenomeAlignment job's data directory, sampleID, BAM path, and
            a flag indicating whether or not the user provided a pre-aligned
            BAM file.

        """
        validation_attributes = {}
        validation_attributes['data_directory'] = self.data_directory_path
        validation_attributes['sample_id'] = sample.sample_id
        if sample.bam_file_path:
            validation_attributes['bam_prealigned'] = True
        else:
            validation_attributes['bam_prealigned'] = False
        validation_attributes['bam_path'] = self.get_genome_bam_path(sample)
        return validation_attributes


    @staticmethod
    def main(cmd_args):
        """
        Entry point into class. Used when script be executed/submitted via the
        command line with the 'align' subcommand.
        """

        sample = eval(cmd_args.sample)
        parameters = json.loads(cmd_args.star_parameters)
        genome_alignment = GenomeAlignmentStep(log_directory_path=cmd_args.log_directory_path,
                                               data_directory_path=cmd_args.data_directory_path,
                                               parameters=parameters)
        genome_alignment.execute(sample=sample,
                                 star_index_directory_path=cmd_args.star_index_directory_path,
                                 star_bin_path=cmd_args.star_bin_path)

    @staticmethod
    def is_output_valid(validation_attributes):
        """
        Check if output of GenomeAlignment for a specific job/execution is
        correctly formed and valid, given a job's data directory, sample id,
        BAM file path, and a flag indicating whether or not the user provided a
        pre-aligned BAM file. If the user provided a pre-aligned BAM file, this
        method assumes that the BAM file is complete if it exists. If this script
        performed the alignment, it will check STAR log files to confirm the BAM
        file is complete.

        Parameters
        ----------
        validation_attributes : dict
            A job's data_directory, sample_id, path to the BAM file, and a flag
            indicating whether or not the user provided a pre-aligned BAM.

        Returns
        -------
        boolean
            True  - GenomeAlignment output files were created and are well formed.
            False - GenomeAlignment output files do not exist or are missing data.
        """
        data_directory = validation_attributes['data_directory']
        sample_id = validation_attributes['sample_id']
        bam_path = validation_attributes['bam_path']
        bam_prealigned = validation_attributes['bam_prealigned']

        valid_output = False

        if os.path.isfile(bam_path):
            if not bam_prealigned:
                alignment_logfile_path = os.path.join(data_directory, f"sample{sample_id}",
                                                      f"{CAMPAREE_CONSTANTS.DEFAULT_STAR_OUTPUT_PREFIX}Log.progress.out")
                if os.path.isfile(alignment_logfile_path):
                    #Read last line in STAR log file
                    line = ""
                    with open(alignment_logfile_path, "r") as alignment_log_file:
                        for line in alignment_log_file:
                            line = line.rstrip()
                        if line == "ALL DONE!":
                            valid_output = True
            else:
                valid_output = True

        return valid_output


class GenomeBamIndexStep(AbstractCampareeStep):

    def __init__(self, log_directory_path, data_directory_path, parameters=None):
        self.log_directory_path = log_directory_path
        self.data_directory_path = data_directory_path

    def validate(self):
        return True

    def execute(self, sample, bam_file_path):
        """
        Build index of a given bam file.

        Parameters
        ----------
        sample : Sample
            Sample associated with BAM file to be indexed.
        bam_file_path : string
            BAM file to be indexed.

        """
        # indexing is necessary for the variant finding step

        if not os.path.isfile(bam_file_path):
            raise CampareeException(f"Cannot index BAM file because it does not"
                                    f" exist: {bam_file_path}")

        print(f"Creating BAM index for sample {sample.sample_name}")
        pysam.index(bam_file_path)
        print(f"Finished creating BAM index for sample {sample.sample_name}")

    def get_commandline_call(self, sample, bam_file_path):
        """
        Prepare command to execute the GenomeIndex from the command line, given
        all of the arugments used to run the execute() function.

        Parameters
        ----------
        sample : Sample
            Sample associated with BAM file to be indexed.
        bam_file_path : string
            BAM file to be indexed.

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.

        """
        #Retrieve path to the genome_alignment.py script.
        genome_align_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        genome_align_path = genome_align_path.rstrip('c')

        command = (f" python {genome_align_path} index"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --sample '{repr(sample)}'"
                   f" --bam_file_path {bam_file_path}")

        return command

    def get_validation_attributes(self, sample, bam_file_path):
        """
        Prepare attributes required by is_output_valid() function to validate
        output generated the BAM index job corresponding to the given bam file.

        Parameters
        ----------
        sample : Sample
            Sample associated with BAM file to be indexed. [Note: this parameter
            is captured just so get_validation_attributes() accepts the same
            arguments as get_commandline_call(). It is not used here.]
        bam_file_path : string
            BAM file to be indexed.

        Returns
        -------
        dict
            Path to the BAM file indexed by this step.

        """
        validation_attributes = {}
        validation_attributes['bam_file_path'] = bam_file_path
        return validation_attributes

    @staticmethod
    def main(cmd_args):
        """
        Entry point into class. Used when script be executed/submitted via the
        command line with the 'index' subcommand.
        """
        sample = eval(cmd_args.sample)
        genome_index = GenomeBamIndexStep(log_directory_path=cmd_args.log_directory_path,
                                          data_directory_path=cmd_args.data_directory_path)
        genome_index.execute(sample=sample,
                             bam_file_path=cmd_args.bam_file_path)

    @staticmethod
    def is_output_valid(validation_attributes):
        """
        Check if output of GenomeBamIndexStep for a specific job/execution is
        valid, given a job's BAM file path.

        Parameters
        ----------
        validation_attributes : dict
            The path to a job's BAM file.

        Returns
        -------
        boolean
            True  - BAM index file was created in same directory as the BAM file.
            False - BAM index file is missing from same directory as the BAM file.

        """

        bam_file_path = validation_attributes['bam_file_path']
        bam_index_file = bam_file_path + ".bai"
        valid_output = False

        if os.path.isfile(bam_index_file):
            valid_output = True

        return valid_output


if __name__ == '__main__':
    """
    Prepare and process command line arguments. The setup belows allows for entry
    into either the GenomeAlignmentStep main() function or the GenomeBamIndexStep
    main() function based on which subcommand is specified at the command line.
    """

    parser = argparse.ArgumentParser(description='Command line wrapper around'
                                                 ' STAR genome alignments and BAM'
                                                 ' file indexing.')

    subparsers = parser.add_subparsers(help="Choose one of the following:",dest="RUN_MODE", metavar="RUN_MODE")
    subparsers.required = True

    #Setup arguments for the alignment subcommand
    genome_alignment_subparser = subparsers.add_parser('align', help="Perform STAR alignment to genome.",
                                                       description="Perform STAR alignment to given genome.")
    genome_alignment_subparser.set_defaults(func=GenomeAlignmentStep.main)
    #Send arguments for this subcommand to the GenomeAlignmentStep's main() method.
    required_named_genome_alignment_subparser = genome_alignment_subparser.add_argument_group('Required named arguments')
    required_named_genome_alignment_subparser.add_argument("--log_directory_path", required=True, help="Directory in which to save logging files.")
    required_named_genome_alignment_subparser.add_argument("--data_directory_path", required=True, help="Directory in which to save output files.")
    required_named_genome_alignment_subparser.add_argument('--star_index_directory_path', required=True, help="Path to directory containing STAR index.")
    required_named_genome_alignment_subparser.add_argument('--sample', required=True, help="Serialized Sample object containing FASTQ file info.")
    required_named_genome_alignment_subparser.add_argument('--star_bin_path', required=True, help="Path to STAR executable binary.")
    required_named_genome_alignment_subparser.add_argument('--star_parameters', required=True, help="Jsonified STAR parameters.")

    #Setup arguments from the index subcommand
    genome_index_subparser = subparsers.add_parser('index', help="Index BAM file aligned to genome.",
                                                   description="Index given BAM file.")
    #Send arguments for this subcommand to the GenomeBamIndexStep's main() method.
    genome_index_subparser.set_defaults(func=GenomeBamIndexStep.main)
    required_named_genome_index_subparser = genome_index_subparser.add_argument_group('Required named arguments')
    required_named_genome_index_subparser.add_argument("--log_directory_path", required=True, help="Directory in which to save logging files.")
    required_named_genome_index_subparser.add_argument("--data_directory_path", required=True, help="Directory in which to save output files.")
    required_named_genome_index_subparser.add_argument('--sample', required=True, help="Serialized Sample object.")
    required_named_genome_index_subparser.add_argument('--bam_file_path', required=True, help="BAM file to be indexed.")

    args = parser.parse_args()
    args.func(args)
