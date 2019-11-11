import argparse
import sys
import os
import collections

from camparee.abstract_camparee_step import AbstractCampareeStep
from camparee.camparee_constants import CAMPAREE_CONSTANTS

class TranscriptGeneQuantificationStep(AbstractCampareeStep):
    """This class takes a kallisto output file and generates transcript- and
    gene-level quantification files, and a file of PSI (Percent Spliced In)
    values for alternative spliceforms.

    """

    OUTPUT_TRANSCRIPT_FILE_NAME = CAMPAREE_CONSTANTS.TXQUANT_OUTPUT_TX_FILENAME
    OUTPUT_GENE_FILE_NAME = CAMPAREE_CONSTANTS.TXQUANT_OUTPUT_GENE_FILENAME
    OUTPUT_PSI_VALUE_FILE_NAME = CAMPAREE_CONSTANTS.TXQUANT_OUTPUT_PSI_FILENAME

    def __init__(self, log_directory_path, data_directory_path):
        """Constructor for TranscriptGeneQuantificationStep object.

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

    def execute(self, sample_id, tx_abundance_file_path, annotation_file_path):
        """Main work-horse function that generates transcript, gene, and PSI
        count files from transcript-level kallisto data.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to the input kallisto file. Used
            to construct output and log paths for this specific execution.
        tx_abundance_file_path : string
            File of transcript abundances created by kallisto. Likely generated
            by KallistoQuantStep.
        annotation_file_path : string
            Input transcript annotation file. Used to map transcript IDs to gene
            IDs. This is generally prepared by the UpdateAnnotationForGenomeStep.

        """

        # Used in create_transcript_gene_map function, so it needs to be an
        # instance variable.
        self.annotation_file_path = annotation_file_path

        # Dictionaries to keep track of length of transcript, number of uniquely
        # mapped reads to transcript, and final count of reads mapped to transcript.
        # This procedure does not map all gene info keys used. Consequently we
        # need to insure that assignments using new keys are initialized to 0.
        transcript_final_count = collections.defaultdict(float)
        psi_value_map = collections.defaultdict(list)
        # Used in create_transcript_gene_map, so it needs to be an instance variable.
        self.transcript_gene_map = collections.defaultdict(str)

        # Prepare paths to output filenames
        transcript_count_filename = os.path.join(self.data_directory_path, f'sample{sample_id}',
                                                 self.OUTPUT_TRANSCRIPT_FILE_NAME)
        gene_count_filename = os.path.join(self.data_directory_path, f'sample{sample_id}',
                                           self.OUTPUT_GENE_FILE_NAME)
        psi_value_filename = os.path.join(self.data_directory_path, f'sample{sample_id}',
                                               self.OUTPUT_PSI_VALUE_FILE_NAME)
        log_file_path = os.path.join(self.log_directory_path, f'sample{sample_id}',
                                     CAMPAREE_CONSTANTS.TXQUANT_LOG_FILENAME)

        with open(log_file_path, "w") as log_file:

            print(f"Generate transcript, gene, and PSI value files from kallisto "
                  f"output from sample{sample_id}.")
            log_file.write("Generate transcript, gene, and PSI value files "
                           f"from kallisto output from sample{sample_id}.\n")
            log_file.write(f"Parameters:\n"
                           f"    kallisto abundance path: {tx_abundance_file_path}\n"
                           f"    annotation file path: {annotation_file_path}\n")

            print("Extracting transcript:gene mappings from annotation file.")
            log_file.write("Extracting transcript:gene mappings from annotation file.\n")
            # Transcript : parent gene dictionary
            self.create_transcript_gene_map()

            print("Extracting transcript abundances from kallisto file.")
            log_file.write("Extracting transcript abundances from kallisto file.\n")
            with open(tx_abundance_file_path, "r") as tx_abundance_file:
                for line in tx_abundance_file:
                    if line.startswith("target_id"):
                        continue
                    line = line.rstrip('\n').split('\t')
                    transcript_id = line[0].split(':')[0]
                    transcript_fpk = float(line[3]) / float(line[2]) * 1000 # est_counts / eff_length * 1000 = FPK
                    transcript_final_count[transcript_id] = transcript_fpk

            print("Writing transcript-level quants to file.")
            log_file.write("Writing transcript-level quants to file.\n")
            # Write the transcript quantification information to transcript quant filename
            with open(transcript_count_filename, 'w') as transcript_count_file:
                transcript_count_file.write('#transcript_id' + '\t' + 'cnt' + '\n')
                for key, value in transcript_final_count.items():
                    transcript_count_file.write(str(key) + '\t' + str(round(value, 3)) + '\n')

            print("Summing transcript-level quants by gene to get gene-level quants.")
            log_file.write("Summing transcript-level quants by gene to get gene-level quants.\n")
            # Add the transcript counts to parent gene to get gene count
            gene_count = collections.defaultdict(float)
            for key, value in transcript_final_count.items():
                gene_count[self.transcript_gene_map[key]] += value

            gene_count = collections.OrderedDict(sorted(gene_count.items()))

            print("Writing gene-level quants to file.")
            log_file.write("Writing gene-level quants to file.\n")
            # Write gene quantification information to gene quant filename
            with open(gene_count_filename, 'w') as gene_count_file:
                gene_count_file.write('#gene_id' + '\t' + 'cnt' + '\n')
                for key, value in gene_count.items():
                    gene_count_file.write(str(key) + '\t' + str(round(value,3)) + '\n')

            print("Deriving PSI values from gene- and transcript-level quants.")
            log_file.write("Deriving PSI values from gene- and transcript-level quants.\n")
            # Create dictionary with key gene_id and values isoforms and their psi values
            for transcript_id in transcript_final_count.keys():
                gene_id = self.transcript_gene_map[transcript_id]
                if gene_count[gene_id]:
                    psi_value_map[gene_id].append(transcript_id + ':' \
                        + str(transcript_final_count[transcript_id]/gene_count[gene_id]))
                else:
                    psi_value_map[gene_id].append(transcript_id + ':' + str(0))

            print("Writing PSI values to file.")
            log_file.write("Writing PSI values to file.\n")
            # Write psi value information for each gene
            with open(psi_value_filename, 'w') as psi_value_file:
                psi_value_file.write('#gene_id' + '\t' + 'isoform_psi_value' + '\n')
                for gene_id in psi_value_map.keys():
                    psi_value_file.write(str(gene_id) + '\t' + ','.join(psi_value_map[gene_id]) + '\n')

            log_file.write("\nALL DONE!\n")

    def create_transcript_gene_map(self):
        """Create dictionary to map transcript id to gene id using geneinfo file

        """
        with open(self.annotation_file_path, 'r') as annotation_file:
            next(annotation_file)
            for line in annotation_file:
                fields = line.strip('\n').split('\t')
                self.transcript_gene_map[fields[7]] = fields[8]

    def get_commandline_call(self, sample_id, tx_abundance_file_path, annotation_file_path):
        """
        Prepare command to execute the TranscriptGeneQuantificationStep from the
        command line, given all of the arugments used to run the execute() function.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to the input kallisto file. Used
            to construct output and log paths for this specific execution.
        tx_abundance_file_path : string
            File of transcript abundances created by kallisto. Likely generated
            by KallistoQuantStep.
        annotation_file_path : string
            Input transcript annotation file. Used to map transcript IDs to gene
            IDs. This is generally prepared by the UpdateAnnotationForGenomeStep.

        Returns
        -------
        string
            Command to execute on the command line. It will perform the same
            operations as a call to execute() with the same parameters.

        """
        #Retrieve path to the transcript_gene_quant.py script.
        txquant_step_path = os.path.realpath(__file__)
        #If the above command returns a string with a "pyc" extension, instead
        #of "py", strip off "c" so it points to this script.
        txquant_step_path = txquant_step_path.rstrip('c')

        command = (f" python {txquant_step_path}"
                   f" --log_directory_path {self.log_directory_path}"
                   f" --data_directory_path {self.data_directory_path}"
                   f" --sample_id {sample_id}"
                   f" --tx_abundance_file_path {tx_abundance_file_path}"
                   f" --annotation_file_path {annotation_file_path}")

        return command

    def get_validation_attributes(self, sample_id, tx_abundance_file_path, annotation_file_path):
        """
        Prepare attributes required by is_output_valid() function to validate
        output generated by the TranscriptGeneQuantificationStep job.

        Parameters
        ----------
        sample_id : string
            Identifier for sample corresponding to the input kallisto file. Used
            to construct output and log paths for this specific execution.
        tx_abundance_file_path : string
            File of transcript abundances created by kallisto. Likely generated
            by KallistoQuantStep. [Note: this parameter is captured just so
            get_validation_attributes() accepts the same arguments as
            get_commandline_call(). It is not used here.]
        annotation_file_path : string
            Input transcript annotation file. Used to map transcript IDs to gene
            IDs. This is generally prepared by the UpdateAnnotationForGenomeStep.
            [Note: this parameter is captured just so get_validation_attributes()
            accepts the same arguments as get_commandline_call(). It is not used
            here.]

        Returns
        -------
        dict
            A TranscriptGeneQuantificationStep job's data_directory, log_directory,
            and corresponding sample ID.
        """
        validation_attributes = {}
        validation_attributes['data_directory'] = self.data_directory_path
        validation_attributes['log_directory'] = self.log_directory_path
        validation_attributes['sample_id'] = sample_id
        return validation_attributes

    @staticmethod
    def is_output_valid(validation_attributes):
        """
        Check if output of TranscriptGeneQuantificationStep for a specific
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
            True  - TranscriptGeneQuantificationStep output files were created
                    and are well formed.
            False - TranscriptGeneQuantificationStep output files do not exist
                    or are missing data.

        """

        data_directory_path = validation_attributes['data_directory']
        log_directory_path = validation_attributes['log_directory']
        sample_id = validation_attributes['sample_id']

        valid_output = False

        # Construct output filenames/paths
        transcript_count_filename = os.path.join(data_directory_path, f'sample{sample_id}',
                                                 TranscriptGeneQuantificationStep.OUTPUT_TRANSCRIPT_FILE_NAME)
        gene_count_filename = os.path.join(data_directory_path, f'sample{sample_id}',
                                           TranscriptGeneQuantificationStep.OUTPUT_GENE_FILE_NAME)
        psi_value_filename = os.path.join(data_directory_path, f'sample{sample_id}',
                                          TranscriptGeneQuantificationStep.OUTPUT_PSI_VALUE_FILE_NAME)
        log_file_path = os.path.join(log_directory_path, f'sample{sample_id}',
                                     CAMPAREE_CONSTANTS.TXQUANT_LOG_FILENAME)

        if os.path.isfile(transcript_count_filename) and \
           os.path.isfile(gene_count_filename) and \
           os.path.isfile(psi_value_filename) and \
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
        """
        Entry point into script. Parses the argument list to obtain all the files
        needed and feeds them to the class constructor. Calls the appropriate
        methods thereafter.

        """

        parser = argparse.ArgumentParser(description='Parse kallisto results and '
                                                     'generate transcript, gene, '
                                                     'and PSI quant files.')
        parser.add_argument('-l', '--log_directory_path', required=True,
                            help="Path to log directory.")
        parser.add_argument('-d', '--data_directory_path', required=True,
                            help='Path to data directory')
        parser.add_argument('--sample_id', required=True,
                            help='Sample ID associated with input kallisto data.')
        parser.add_argument('-k', '--tx_abundance_file_path', required=True,
                            help='Transcript-level abundance file generated by kallisto')
        parser.add_argument('-a', '--annotation_file_path', required=True,
                            help='Annotation file mapping transcript IDs to gene IDs')

        args = parser.parse_args()

        transcript_gene_quant = TranscriptGeneQuantificationStep(log_directory_path=args.log_directory_path,
                                                                 data_directory_path=args.data_directory_path)
        transcript_gene_quant.execute(sample_id=args.sample_id,
                                      tx_abundance_file_path=args.tx_abundance_file_path,
                                      annotation_file_path=args.annotation_file_path)

if __name__ == "__main__":
    sys.exit(TranscriptGeneQuantificationStep.main())
