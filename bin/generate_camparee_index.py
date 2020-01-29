#!/usr/bin/env python
"""
This module generates the reference genome sequence and annotation files in the
format required by CAMPAREE. As input, it requires the reference sequence in
FASTA format (one entry for each contig/chromosome), the annotation file in GTF
format, and the name of the species/model corresponding to the input files (e.g.
Mus_musculus.Ensembl_v97.GRCm38). The output files will be stored in the
resources/ directory of this CAMPAREE install, under the species/model name
given as part of the input.
"""
import argparse
import sys
import os

from camparee.camparee_constants import CAMPAREE_CONSTANTS
from camparee.camparee_utils import CampareeUtils

class GenerateCampareeIndex():

    def __init__(self, species_model, genome_fasta_filename, annotation_gtf_filename, output_directory_path):
        """Short summary.

        Parameters
        ----------
        species_model : string
            Description of parameter `species_model`.
        genome_fasta_filename : string
            Description of parameter `genome_fasta_filename`.
        annotation_gtf_filename : string
            Description of parameter `annotation_gtf_filename`.
        output_directory_path: string
            Description of parameter `output_directory_path`.

        """
        self.species_model = species_model
        self.genome_fasta_filename = genome_fasta_filename
        self.annotation_gtf_filename = annotation_gtf_filename

        # Test species_model for no spaces or special characters

        self.output_directory = ""
        if not output_directory_path:
            self.output_directory = os.path.join(CAMPAREE_CONSTANTS.CAMPAREE_ROOT_DIR,
                                                 "resources", self.species_model)
        else:
            # Check if upper-level output_directory_path exists
            if os.path.isdir(output_directory_path):
                self.output_directory = os.path.join(output_directory_path, self.species_model)
            else:
                raise Exception(f"The given output directory does not exist:\n"
                                f"{output_directory_path}")

        # Create output directory if it does not already exist.
        if not os.path.isdir(self.output_directory):
            os.mkdir(self.output_directory)

        # Test for existence of input files (if they were provided).
        if genome_fasta_filename and not os.path.isfile(self.genome_fasta_filename):
            sys.exit(f"ERROR: The given genome FASTA file does not exist:\n"
                     f"    {self.genome_fasta_filename}\n")
        if annotation_gtf_filename and not os.path.isfile(self.annotation_gtf_filename):
            sys.exit(f"ERROR: The given annotation GTF file does not exist:\n"
                     f"    {self.annotation_gtf_filename}\n")

    # Method to format FASTA file
    def format_fasta_file_for_camparee(self):

        output_fasta_file_path = os.path.join(self.output_directory,
                                              self.species_model + ".oneline_seqs.fa")
        # Check if output file already exists. Prevents accidental overwriting
        # of index files.
        if os.path.isfile(output_fasta_file_path):
            sys.exit(f"ERROR: Name of output FASTA file formatted for CAMPAREE "
                     f"already exists:\n    {output_fasta_file_path}\n")

        CampareeUtils.create_oneline_seq_fasta(self.genome_fasta_filename, output_fasta_file_path)

    # Method to format GTF file
    def format_gtf_file_for_camparee(self):

        output_annot_file_path = os.path.join(self.output_directory,
                                              self.species_model + ".annotation.txt")

        # Check if output file already exists. Prevents accidental overwriting
        # of index files.
        if os.path.isfile(output_annot_file_path):
            sys.exit(f"ERROR: Output annotation file formatted for CAMPAREE "
                     f"already exists:\n    {output_annot_file_path}\n")

        CampareeUtils.convert_gtf_to_annot_file_format(input_gtf_filename=self.annotation_gtf_filename,
                                                       output_annot_filename=output_annot_file_path)

    # Maybe method to generate STAR index?

    @staticmethod
    def main():
        """Entry point into script.

        Parses arguments, gathers input filenames, and calls methods that perform
        the actual operation.

        """
        parser = argparse.ArgumentParser(description='Format genome FASTA and/or annotation'
                                                     ' GTF files for use with CAMPAREE.')
        parser.add_argument('-n', '--species_model_name', required=True,
                            help="Name of the species and genome build of the input files. "
                                 "Used to name output directory and files, so it must be "
                                 "compatible with file naming conventions.")
        parser.add_argument('-g', '--genome_fasta_file', required=False, default="",
                            help="Path to genome sequence in FASTA format.")
        parser.add_argument('-a', '--annotation_gtf_file', required=False, default="",
                            help="Path to annotation file in GTF format.")
        parser.add_argument('-o', '--output_directory_path', required=False, default="",
                            help="Path to output the formatted files. If not specified, "
                                 "files are saved to resources/ subdirectory within "
                                 "CAMPAREE install directory.")
        args = parser.parse_args()

        if not (args.genome_fasta_file or args.annotation_gtf_file):
            parser.error('No operation performed. Please enter a FASTA and/or a GTF file '
                         'to format for CAMPAREE.')

        camparee_index_generator = GenerateCampareeIndex(species_model=args.species_model_name,
                                                         genome_fasta_filename=args.genome_fasta_file,
                                                         annotation_gtf_filename=args.annotation_gtf_file,
                                                         output_directory_path=args.output_directory_path)

        if args.genome_fasta_file:
            camparee_index_generator.format_fasta_file_for_camparee()
        if args.annotation_gtf_file:
            camparee_index_generator.format_gtf_file_for_camparee()




if __name__ == '__main__':
    sys.exit(GenerateCampareeIndex.main())
