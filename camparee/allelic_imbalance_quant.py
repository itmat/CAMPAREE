import argparse
import re
import sys
import os
import collections

from pysam import AlignmentFile

from camparee.abstract_camparee_step import AbstractCampareeStep
from camparee.camparee_constants import CAMPAREE_CONSTANTS

# TODO: Go back through and optimize this code to use fewer class variables
#       (could pass necessary info as arguments to helper functions).

class AllelicImbalanceQuantificationStep(AbstractCampareeStep):
    """This class contains scripts to output quantification of allelic imbalance.

    It requires
     (i) an input file source for gene info
    (ii) Root of the aligned filenames (alignment to transcriptome of each parent
         with suffixes '_1','_2'.)

    There is one output file with quantification information on the allelic
    imbalance of genes. Fields in this file: chromosome, strand, start, end,
    exon count, exon starts, exon ends, gene name.

    """

    OUTPUT_ALLELIC_IMBALANCE_FILE_NAME = CAMPAREE_CONSTANTS.ALLELIC_IMBALANCE_OUTPUT_FILENAME

    def __init__(self, log_directory_path, data_directory_path, parameters=None):
        """Constructor for AllelicImbalanceQuantificationStep object.

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

    def create_transcript_gene_map(self):
        """
        Create dictionary to map transcript id to gene id using geneinfo file
        Map '*' to '*' to account for unmapped reads in align_file
        Create entries with suffix '_1' and '_2' for each transcript

        """
        self.transcript_gene_map['*'] = '*'

        with open(self.geneinfo_filename_1, 'r') as geneinfo_file:
            next(geneinfo_file)
            for line in geneinfo_file:
                fields = line.strip('\n').split('\t')
                self.transcript_gene_map[fields[7]] = fields[8]

    def reads_to_ignore(self):
        reads_to_ignore = []
        bamfile = AlignmentFile(self.genome_alignment_file, "rb")
        num_hits_pattern = re.compile('(NH:i:)(\d+)')

        for read in bamfile.fetch(until_eof=True):
            num_hits = dict(read.tags)['NH']
            if num_hits > 1:
                reads_to_ignore.append(read.query_name)

        return reads_to_ignore

    def read_info(self, in_align_filename):
        """
        Create dictionary which maps a read id in SAM file to a dictionary with two keys 'transcript_id' and 'NM'.
        The value associated with 'transcript_id' is a list of all transcripts the read aligned to.
        The value associated with 'NM' is the corresponding edit distance information for each alignment.
        For non-mappers the transcript_id is '*' and edit distance is 100 (Make it read length).

        """
        read_info_map = collections.defaultdict(dict)

        # The NM tag in the SAM file tells us the edit distance for the alignment.
        # This pattern extracts that number.
        num_mismatches_pattern = re.compile('(NM:i:)(\d+)')

        with open(in_align_filename, 'r') as infile:
            for line in infile:
                if line.startswith('@'):
                    continue

                # read forward and reverse read
                forward = line
                reverse = next(infile)

                # Parse the fields for the forward read into an array
                fwd_fields = forward.rstrip('\n').split('\t')

                # Parse the fields for the reverse read into an array
                rev_fields = reverse.rstrip('\n').split('\t')

                fwd_transcript_id = fwd_fields[2].split(':')[0]
                rev_transcript_id = rev_fields[2].split(':')[0]

                # This means both forward and reverse reads are non-mappers
                # So store 'transcript_id' as '*' and 'NM' as 2*read_length
                if fwd_transcript_id == '*' and rev_transcript_id == '*':
                    read_info_map[fwd_fields[0]]['transcript_id'] =  '*'
                    read_info_map[fwd_fields[0]]['NM'] =  200
                    continue
                # Get transcript_id for mapped reads
                elif fwd_transcript_id == rev_transcript_id:
                    transcript_id = fwd_transcript_id
                else:
                    transcript_id = (fwd_transcript_id + rev_transcript_id).replace('*','')

                # This probably means the transcript was not in our master list of all transcript models
                #  (the geneinfo filename).  So we skip it.  Really this should not happen
                #  but just in case.
                if not self.transcript_gene_map.get(transcript_id):
                    continue

                # Obtain the edit distance information for the forward read
                fwd_NM_match = re.search(num_mismatches_pattern, forward)
                rev_NM_match = re.search(num_mismatches_pattern, reverse)
                if fwd_NM_match and rev_NM_match:
                    fwd_NM_count = int(fwd_NM_match.group(2))
                    rev_NM_count = int(rev_NM_match.group(2))
                    NM_count = fwd_NM_count + rev_NM_count
                elif not (fwd_NM_match and rev_NM_match):
                    NM_count = 200
                else:
                    NM_count = 100

                # Update read_info dictionary with transcript_id and corresponding edit distance
                read_info_map[fwd_fields[0]]['transcript_id'] = transcript_id
                read_info_map[fwd_fields[0]]['NM'] = NM_count

        return read_info_map

    def execute(self, sample_id, genome_alignment_file_path, parent1_annot_file_path,
                parent2_annot_file_path, parent1_tx_align_file_path, parent2_tx_align_file_path):
        """This is the main method which quantifies allelic imbalance for all
        genes in the annotation based on the aligned files for parents 1 and 2.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to the input genome and transcriptome
            alignment files. Used to construct output and log paths for this specific
            execution.
        genome_alignment_file_path : string
            Input BAM file of reads aligned to the original reference genome.
            This is used to identify multimappers so they are excluded from the
            allelic imbalance quantification. This is generally prepared by
            GenomeAlignmentStep, or provided by the user.
        parent1_annot_file_path : string
            Input transcript annotation file for parent 1. This is generally
            prepared by UpdateAnnotationForGenomeStep.
        parent2_annot_file_path : string
            Input transcript annotation file for parent 2. This is generally
            prepared by UpdateAnnotationForGenomeStep.
        parent1_tx_align_file_path : string
            Input SAM file of reads aligned to the variant genome from parent 1.
            This is generally prepared by Bowtie2AlignStep.
        parent2_tx_align_file_path : string
            Input SAM file of reads aligned to the variant genome from parent 2.
            This is generally prepared by Bowtie2AlignStep.

        """
        self.genome_alignment_file = genome_alignment_file_path
        self.geneinfo_filename_1 = parent1_annot_file_path
        self.geneinfo_filename_2 = parent2_annot_file_path
        self.align_filename_1 = parent1_tx_align_file_path
        self.align_filename_2 = parent2_tx_align_file_path

        log_file_path = os.path.join(self.log_directory_path, f'sample{sample_id}',
                                     CAMPAREE_CONSTANTS.ALLELIC_IMBALANCE_LOG_FILENAME)

        # Create allelic imbalance distribution file and ensure that it doesn't
        # currently exist.
        self.allele_imbalance_dist_filename = os.path.join(self.data_directory_path, f'sample{sample_id}',
                                                           AllelicImbalanceQuantificationStep.OUTPUT_ALLELIC_IMBALANCE_FILE_NAME)
        try:
            os.remove(self.allele_imbalance_dist_filename)
        except OSError:
            pass

        # Dictionaries to map transcripts to genes and keep track of final count of reads mapped to a gene
        # This procedure does not map all gene info keys used.  Consequently we need to insure that
        # assignments using new keys are initialized to 0
        self.transcript_gene_map = collections.defaultdict(str)
        self.gene_final_count = collections.defaultdict(lambda: collections.defaultdict(int))
        self.exclusive_genes = []

        with open(log_file_path, "w") as log_file:

            print(f"Quantify allelic imbalance for reads from sample{sample_id}.")
            log_file.write(f"Quantify allelic imbalance for reads from sample{sample_id}.\n")
            log_file.write(f"Parameters:\n"
                           f"    Reference genome align path:       {self.genome_alignment_file}\n"
                           f"    Parent 1 annotation path:          {self.geneinfo_filename_1}\n"
                           f"    Parent 2 annotation path:          {self.geneinfo_filename_2}\n"
                           f"    Parent 1 transcriptome align path: {self.align_filename_1}\n"
                           f"    Parent 2 transcriptome align path: {self.align_filename_2}\n")

            print("Mapping transcript IDs to gene IDs from Parent 1 annotation file.")
            log_file.write("Mapping transcript IDs to gene IDs from Parent 1 annotation file.\n")
            self.create_transcript_gene_map()

            # Create read info dictionaries for each parent
            print("Extracting read-transcript mappings for parent 1 from"
                  " transcriptome alignment file.")
            log_file.write("Extracting read-transcript mappings for parent 1"
                           " from transcriptome alignment file.")
            read_info_1 = self.read_info(self.align_filename_1)
            print("Extracting read-transcript mappings for parent 2 from"
                  " transcriptome alignment file.")
            log_file.write("Extracting read-transcript mappings for parent 2"
                           " from transcriptome alignment file.")
            read_info_2 = self.read_info(self.align_filename_2)

            print("Identifying multimappers from genome alignments.")
            log_file.write("Identifying multimappers from genome alignments.\n")
            reads_to_ignore = self.reads_to_ignore()

            print("Excluding multimappers from further use.")
            log_file.write("Excluding multimappers from further use.\n")
            read_ids_1 = set(read_info_1.keys()).difference(reads_to_ignore)
            read_ids_2 = set(read_info_2.keys()).difference(reads_to_ignore)

            read_ids = read_ids_1.intersection(read_ids_2)
            read_ids_1_u = read_ids_1.difference(read_ids)
            read_ids_2_u = read_ids_2.difference(read_ids)

            print("Quantifying reads aligned only to parental genome 1.")
            log_file.write("Quantifying reads aligned only to parental genome 1.\n")
            for read in read_ids_1_u:
                transcript = read_info_1[read]['transcript_id']
                gene = self.transcript_gene_map[transcript]
                self.gene_final_count[gene]['1'] += 1

            print("Quantifying reads aligned only to parental genome 2.")
            log_file.write("Quantifying reads aligned only to parental genome 2.\n")
            for read in read_ids_2_u:
                transcript = read_info_2[read]['transcript_id']
                gene = self.transcript_gene_map[transcript]
                self.gene_final_count[gene]['2'] += 1

            print("Quantifying reads aligned to both parental genomes.")
            log_file.write("Quantifying reads aligned to both parental genomes.\n")
            for read in read_ids:
                # Transcripts to which the read mapped for each parent
                transcript_1 = read_info_1[read]['transcript_id']
                transcript_2 = read_info_2[read]['transcript_id']

                # The read did not map to any transcript in either parent
                if transcript_1 == '*' and transcript_2 == '*':
                    continue
                # The read mapped to atleast one transcript in each parent
                elif transcript_1 != '*' and transcript_2 != '*':
                    # Get the genes in parent 1 to which the read mapped
                    gene_1 = self.transcript_gene_map[transcript_1]
                    NM_count_1 = read_info_1[read]['NM']

                    # Get the genes in parent 2 to which the read mapped
                    gene_2 = self.transcript_gene_map[transcript_2]
                    NM_count_2 = read_info_2[read]['NM']

                    # Amongst the genes to which the read mapped,
                    # there is exactly one gene in common between parent 1 and 2.
                    if gene_1 == gene_2:
                        # Minimum edit distance for the mapping to the gene is the same in
                        # parent 1 and parent 2. So increment counts of both alleles of the genes by 0.5
                        if NM_count_1 == NM_count_2:
                            self.gene_final_count[gene_1]['1'] += 0.5
                            self.gene_final_count[gene_1]['2'] += 0.5
                        # Minimum edit distance for the mapping to the gene is less in parent 1.
                        # So increment count of allele of gene corresponding to parent 1.
                        elif NM_count_1 < NM_count_2:
                            self.gene_final_count[gene_1]['1'] += 1
                        # Minimum edit distance for the mapping to the gene is less in parent 2.
                        # So increment count of allele of gene corresponding to parent 2.
                        else:
                            self.gene_final_count[gene_1]['2'] += 1
                # The read is a non-mapper for the parent 1 transcriptome
                elif transcript_1 == '*':
                    # Get the genes in parent 2 to which the read mapped
                    gene_2 = self.transcript_gene_map[transcript_2]
                    self.gene_final_count[gene_2]['2'] += 1

                # The read is a non-mapper for the parent 2 transcriptome
                elif transcript_2 == '*':
                    gene_1 = self.transcript_gene_map[transcript_1]
                    self.gene_final_count[gene_1]['1'] += 1

            print("Writing file of allelic imbalance quantification results.")
            log_file.write("Writing file of allelic imbalance quantification results.\n")
            self.make_allele_imbalance_dist_file()

            log_file.write("\nALL DONE!\n")

    def make_allele_imbalance_dist_file(self):
        genelist_1 = []
        with open(self.geneinfo_filename_1, 'r') as geneinfo_file_1:
            for line in geneinfo_file_1:
                if line.startswith('#'):
                    continue
                fields = line.strip('\n').split('\t')
                genelist_1.append(fields[8])

        genelist_2 = []
        with open(self.geneinfo_filename_2, 'r') as geneinfo_file_2:
            for line in geneinfo_file_2:
                if line.startswith('#'):
                    continue
                fields = line.strip('\n').split('\t')
                genelist_2.append(fields[8])

        exclusive_genes = list(set(genelist_1).difference(set(genelist_2)))

        # Write the allelic imbalance quantification information to allele imbalance dist filename
        with open(self.allele_imbalance_dist_filename, 'w') as allele_imbalance_dist_file:
            allele_imbalance_dist_file.write('#gene_id' + '\t' + '_1' + '\t' + '_2' + '\n')

            #for key, value in list(self.gene_final_count.items()):
            for gene_id in sorted(set(self.transcript_gene_map.values())):
                if gene_id in exclusive_genes:
                    allele_imbalance_dist_file.write(str(gene_id) + '\t' + str(1.0) + '\t' + str(0.0) + '\n')
                    continue

                if gene_id == "*":
                    continue

                read_count_1 = self.gene_final_count[gene_id]['1']
                read_count_2 = self.gene_final_count[gene_id]['2']
                gene_read_count = read_count_1 + read_count_2

                if gene_read_count == 0:
                    allele_imbalance_dist_file.write(str(gene_id) + '\t' + str(0.5) + '\t' + str(0.5) + '\n')
                else:
                    allele_imbalance_dist_file.write(str(gene_id) + '\t' + str(read_count_1/gene_read_count) + '\t' +\
                        str(read_count_2/gene_read_count) + '\n')

    def get_commandline_call(self, sample_id, genome_alignment_file_path,
                             parent1_annot_file_path, parent2_annot_file_path,
                             parent1_tx_align_file_path, parent2_tx_align_file_path):
        """Prepare command to execute the AllelicImbalanceQuantificationStep from
        the command line, given all of the arugments used to run the execute()
        function.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to the input genome and transcriptome
            alignment files. Used to construct output and log paths for this specific
            execution.
        genome_alignment_file_path : string
            Input BAM file of reads aligned to the original reference genome.
            This is used to identify multimappers so they are excluded from the
            allelic imbalance quantification. This is generally prepared by
            GenomeAlignmentStep, or provided by the user.
        parent1_annot_file_path : string
            Input transcript annotation file for parent 1. This is generally
            prepared by UpdateAnnotationForGenomeStep.
        parent2_annot_file_path : string
            Input transcript annotation file for parent 2. This is generally
            prepared by UpdateAnnotationForGenomeStep.
        parent1_tx_align_file_path : string
            Input SAM file of reads aligned to the variant genome from parent 1.
            This is generally prepared by Bowtie2AlignStep.
        parent2_tx_align_file_path : string
            Input SAM file of reads aligned to the variant genome from parent 2.
            This is generally prepared by Bowtie2AlignStep.

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.

        """
        #Retrieve path to the allelic_imbalance_quant.py script.
        allelic_imbalance_step_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        allelic_imbalance_step_path = allelic_imbalance_step_path.rstrip('c')

        command = (f" python {allelic_imbalance_step_path}"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --sample_id {sample_id}"
                   f" --genome_alignment_path {genome_alignment_file_path}"
                   f" --parent1_annot_path {parent1_annot_file_path}"
                   f" --parent2_annot_path {parent2_annot_file_path}"
                   f" --parent1_tx_align_path {parent1_tx_align_file_path}"
                   f" --parent2_tx_align_path {parent2_tx_align_file_path}")

        return command

    def get_validation_attributes(self, sample_id, genome_alignment_file_path,
                             parent1_annot_file_path, parent2_annot_file_path,
                             parent1_tx_align_file_path, parent2_tx_align_file_path):
        """Prepare attributes required by is_output_valid() function to validate
        output generated by the AllelicImbalanceQuantificationStep job.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to the input genome and transcriptome
            alignment files. Used to construct output and log paths for this specific
            execution.
        genome_alignment_file_path : string
            Input BAM file of reads aligned to the original reference genome.
            This is used to identify multimappers so they are excluded from the
            allelic imbalance quantification. This is generally prepared by
            GenomeAlignmentStep, or provided by the user. [Note: this parameter
            is captured just so get_validation_attributes() accepts the same
            arguments as get_commandline_call(). It is not used here.]
        parent1_annot_file_path : string
            Input transcript annotation file for parent 1. This is generally
            prepared by UpdateAnnotationForGenomeStep. [Note: this parameter
            is captured just so get_validation_attributes() accepts the same
            arguments as get_commandline_call(). It is not used here.]
        parent2_annot_file_path : string
            Input transcript annotation file for parent 2. This is generally
            prepared by UpdateAnnotationForGenomeStep. [Note: this parameter
            is captured just so get_validation_attributes() accepts the same
            arguments as get_commandline_call(). It is not used here.]
        parent1_tx_align_file_path : string
            Input SAM file of reads aligned to the variant genome from parent 1.
            This is generally prepared by Bowtie2AlignStep. [Note: this parameter
            is captured just so get_validation_attributes() accepts the same
            arguments as get_commandline_call(). It is not used here.]
        parent2_tx_align_file_path : string
            Input SAM file of reads aligned to the variant genome from parent 2.
            This is generally prepared by Bowtie2AlignStep. [Note: this parameter
            is captured just so get_validation_attributes() accepts the same
            arguments as get_commandline_call(). It is not used here.]

        Returns
        -------
        dict
            A AllelicImbalanceQuantificationStep job's data_directory, log_directory,
            and corresponding sample ID.
        """
        validation_attributes = {}
        validation_attributes['data_directory'] = self.data_directory_path
        validation_attributes['log_directory'] = self.log_directory_path
        validation_attributes['sample_id'] = sample_id
        return validation_attributes

    @staticmethod
    def is_output_valid(validation_attributes):
        """Check if output of AllelicImbalanceQuantificationStep for a specific
        job/execution is correctly formed and valid, given a job's data directory,
        log directory, and sample ID. Prepare these attributes for a given job
        using the get_validation_attributes() method.

        Parameters
        ----------
        validation_attributes : dict
            A job's data_directory, log_directory, and corresponding sample_id.

        Returns
        -------
        boolean
            True  - AllelicImbalanceQuantificationStep output files were created
                    and are well formed.
            False - AllelicImbalanceQuantificationStep output files do not exist
                    or are missing data.

        """

        data_directory_path = validation_attributes['data_directory']
        log_directory_path = validation_attributes['log_directory']
        sample_id = validation_attributes['sample_id']

        valid_output = False

        # Construct output filenames/paths
        allele_imbalance_dist_filename = os.path.join(data_directory_path, f'sample{sample_id}',
                                                      AllelicImbalanceQuantificationStep.OUTPUT_ALLELIC_IMBALANCE_FILE_NAME)
        log_file_path = os.path.join(log_directory_path, f'sample{sample_id}',
                                     CAMPAREE_CONSTANTS.ALLELIC_IMBALANCE_LOG_FILENAME)

        # TODO: Report reason why out validation failed.

        if os.path.isfile(allele_imbalance_dist_filename) and \
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
        """Entry point into script. Parses the argument list to obtain all the
        files needed and feeds them to the class constructor. Calls the appropriate
        methods thereafter.

        """

        parser = argparse.ArgumentParser(description='Generate allelic imbalance quantifications')
        parser.add_argument('-l', '--log_directory_path', required=True,
                            help="Path to log directory.")
        parser.add_argument('-d', '--data_directory_path', required=True,
                            help='Path to data directory')
        parser.add_argument('--sample_id', required=True,
                            help='Sample ID associated with input data.')
        parser.add_argument('--genome_alignment_path', required=True,
                            help='BAM file of reads aligned to reference.')
        parser.add_argument('--parent1_annot_path', required=True,
                            help='Annotation file from genome for parent 1.')
        parser.add_argument('--parent2_annot_path', required=True,
                            help='Annotation file from genome for parent 2.')
        parser.add_argument('--parent1_tx_align_path', required=True,
                            help='SAM file of reads aligned to parent 1 transcriptome.')
        parser.add_argument('--parent2_tx_align_path', required=True,
                            help='SAM file of reads aligned to parent 2 transcriptome.')

        args = parser.parse_args()

        transcript_gene_quant = AllelicImbalanceQuantificationStep(log_directory_path=args.log_directory_path,
                                                                   data_directory_path=args.data_directory_path)
        transcript_gene_quant.execute(sample_id=args.sample_id,
                                      genome_alignment_file_path=args.genome_alignment_path,
                                      parent1_annot_file_path=args.parent1_annot_path,
                                      parent2_annot_file_path=args.parent2_annot_path,
                                      parent1_tx_align_file_path=args.parent1_tx_align_path,
                                      parent2_tx_align_file_path=args.parent2_tx_align_path)

if __name__ == "__main__":
    sys.exit(AllelicImbalanceQuantificationStep.main())
