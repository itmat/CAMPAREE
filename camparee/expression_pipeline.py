import sys
import os
import importlib
import inspect
import numpy

from beers_utils.constants import CONSTANTS,SUPPORTED_SCHEDULER_MODES,MAX_SEED
from camparee.camparee_constants import CAMPAREE_CONSTANTS
from beers_utils.job_monitor import JobMonitor
from camparee.camparee_utils import CampareeUtils, CampareeException

class ExpressionPipeline:
    """
    This class represents a pipeline of steps that take user supplied fastq files through alignment, variants
    finding, parental genome construction, annotation, quantification and generation of transcripts and finally the
    generation of packets of molecules that may be used to simulate RNA sequencing.
    """

    THIRD_PARTY_SOFTWARE_DIR_PATH = os.path.join(CAMPAREE_CONSTANTS.CAMPAREE_ROOT_DIR, "third_party_software")
    REQUIRED_RESOURCE_MAPPINGS = ['species_model', 'star_genome_index_directory_name',
                                  'reference_genome_filename', 'annotation_filename', 'chr_ploidy_filename']
    REQUIRED_OUTPUT_MAPPINGS = ['directory_path', 'type', 'override_sample_molecule_count', 'default_molecule_count']

    def __init__(self, configuration, scheduler_mode, output_directory_path, input_samples):
        self.scheduler_mode = scheduler_mode
        self.scheduler_default_params = {'default_num_processors': None,
                                         'default_memory_in_mb': None,
                                         'default_submission_args': None}
        self.samples = input_samples
        self.output_directory_path = output_directory_path
        log_directory_path = os.path.join(output_directory_path, CONSTANTS.LOG_DIRECTORY_NAME)
        data_directory_path = os.path.join(output_directory_path, CONSTANTS.DATA_DIRECTORY_NAME)
        self.data_directory_path = data_directory_path
        self.log_directory_path = log_directory_path
        self.create_intermediate_data_subdirectories(data_directory_path, log_directory_path)
        self.log_file_path = os.path.join(log_directory_path, "expression_pipeline.log")
        self.optional_inputs = {}
        self.sample_optional_inputs = {}
        self.steps = {}
        #Track pathes of scripts for each step. This is needed when running the
        #steps from the command line, as we do when submitting to lsf.
        self.__step_paths = {}
        # Individual steps can provide scheduler parameters that override
        # the default values.
        self.__step_scheduler_param_overrides = {}

        # Validate the resources and set file and directory paths as needed.
        if not self.validate_and_set_resources(configuration['resources']):
            raise CampareeValidationException("The resources data is not completely valid.  "
                                              "Consult the standard error file for details.")

        # Collect the data from the ref genome and chromosome ploidy files
        self.reference_genome = CampareeUtils.create_genome(self.reference_genome_file_path)
        self.chr_ploidy_data = CampareeUtils.create_chr_ploidy_data(self.chr_ploidy_file_path)

        # Set 3rd party software paths
        self.set_third_party_software()

        # Validate output data and set
        if not self.validate_and_set_output_data(configuration['output']):
            raise CampareeValidationException("The output data is not completely valid.  "
                                              "Consult the standard error file for details.")

        # Validate and instantiate optional input files
        if not self.validate_and_set_optional_inputs(configuration['input']['optional_inputs']):
            raise CampareeValidationException("There was a problem with the optional inputs. "
                                              "Consult the standard error file for details.")

        # Validate and instantiate sample-specific optional input files
        if not self.validate_and_set_sample_optional_inputs(configuration['input']['data']):
            raise CampareeValidationException("There was a problem with the per sample "
                                              "optional inputs. Consult the standard error "
                                              "file for details.")

        # Load default scheduler parameters (if provided)
        print(f"Running CAMPAREE using the {self.scheduler_mode} job scheduler.",
              file=sys.stderr)
        if configuration['setup'].get('default_scheduler_parameters'):
            self.scheduler_default_params['default_num_processors'] = configuration['setup']['default_scheduler_parameters'].get('default_num_processors', None)
            self.scheduler_default_params['default_memory_in_mb'] = configuration['setup']['default_scheduler_parameters'].get('default_memory_in_mb', None)
            self.scheduler_default_params['default_submission_args'] = configuration['setup']['default_scheduler_parameters'].get('default_submission_args', None)
            print("With the following default scheduler parameters:\n",
                  '\n'.join({f"\t-{key} : {value}" for key, value in self.scheduler_default_params.items()}),
                  file=sys.stderr)

        # Validate job resubmission limit (if provided)
        # TODO: Update BEERS_UTILS to specify the max_resub_limit as a constant,
        #       which code here can reference for error output rather than
        #       using a separate, CAMPAREE-specific default value.
        self.max_resub_limit = configuration['setup'].get('job_resub_limit', 3)
        if not isinstance(self.max_resub_limit, int) or self.max_resub_limit < 0:
            print(f"The given job resubmission limit (job_resub_limit={self.max_resub_limit})",
                  "is invalid. It must be an integer value >= 0.",
                  file=sys.stderr)
            raise CampareeValidationException("There was a problem with the optional job "
                                              "resubmission limit (job_resub_limit). "
                                              "Consult the standard error file for details.")
        print(f"And a maximum job resubmission limit of {self.max_resub_limit}.")

        self.expression_pipeline_monitor = JobMonitor(output_directory_path=self.output_directory_path,
                                                      scheduler_name=self.scheduler_mode,
                                                      default_num_processors=self.scheduler_default_params['default_num_processors'],
                                                      default_memory_in_mb=self.scheduler_default_params['default_memory_in_mb'],
                                                      max_resub_limit=self.max_resub_limit)

        # Load instances of each pipeline step into the dictionyar of pipleine
        # steps tracked by the job monitor.
        for step, props in configuration['steps'].items():
            module_name, step_name = step.rsplit(".")
            parameters = props["parameters"] if props and "parameters" in props else None
            scheduler_parameters = props["scheduler_parameters"] if props and "scheduler_parameters" in props else None
            module = importlib.import_module(f'.{module_name}', package="camparee")
            step_class = getattr(module, step_name)
            self.steps[step_name] = step_class(log_directory_path, data_directory_path, parameters)
            self.__step_paths[step_name] = inspect.getfile(module)
            self.__step_scheduler_param_overrides[step_name] = scheduler_parameters
            self.expression_pipeline_monitor.add_pipeline_step(step_name=step_name,
                                                               step_class=step_class)

        # The MoleculeMaker step is configured a little differently from the other
        # steps. Most of the steps are instantiated and configured based on entries
        # in the "steps:" section of the config file. MoleculeMaker does not have
        # an entry in the steps section because it is configured from the "output:"
        # section of the config file.
        module_name = "molecule_maker"
        step_name = "MoleculeMakerStep"
        parameters = None
        scheduler_parameters = None
        module = importlib.import_module(f'.{module_name}', package="camparee")
        step_class = getattr(module, step_name)
        self.steps[step_name] = step_class(log_directory_path, data_directory_path, parameters)
        self.__step_paths[step_name] = inspect.getfile(module)
        self.__step_scheduler_param_overrides[step_name] = scheduler_parameters
        self.expression_pipeline_monitor.add_pipeline_step(step_name, step_class)

    def create_intermediate_data_subdirectories(self, data_directory_path, log_directory_path):
        for sample in self.samples:
            os.makedirs(os.path.join(data_directory_path, f'sample{sample.sample_id}'), mode=0o0755, exist_ok=True)
            os.makedirs(os.path.join(log_directory_path, f'sample{sample.sample_id}'), mode=0o0755, exist_ok=True)

    def validate_and_set_output_data(self, output):
        """
        Helper method to validate and set output data.
        :param output: The output dictionary extracted from the configuration file.
        :return: True for valid output data and False otherwise
        """

        valid = True

        # Insure that all required resources keys are in place.  No point in continuing until this
        # problem is resolved.
        missing_output_keys = [item
                                 for item in ExpressionPipeline.REQUIRED_OUTPUT_MAPPINGS
                                 if item not in output]
        if missing_output_keys:
            print(f"The following required mappings were not found under 'outputs': "
                  f"{(',').join(missing_output_keys)}", file=sys.stderr)
            return False

        # Insure type mapping exists
        if "type" not in output:
            print(f"The required mapping 'type' was not found under 'output.", file=sys.stderr)
            valid = False
        else:
            self.output_type = output["type"]

        # Insure default_molecule_count exists and is an int
        # TODO: is this redundant given the check performed with missing_output_keys above?
        if "default_molecule_count" not in output:
            print(f"The required mapping 'default_molecule_count' was not found under 'output.", file=sys.stderr)
            valid = False
        else:
            self.default_molecule_count = output["default_molecule_count"]
            if not isinstance(self.default_molecule_count, int) and \
               self.default_molecule_count < 0:
                print(f"The 'default_molecule_count' must be a positive integer - not "
                      f"{output['default_molecule_count']}",
                      file=sys.stderr)
                valid = False

        #Check to make sure override_sample_molecule_count is a boolean
        self.override_sample_molecule_count = output["override_sample_molecule_count"]
        if not isinstance(self.override_sample_molecule_count, bool):
            print(f"The 'override_sample_molecule_count' must be a boolean (True/False) "
                  f"- not {output['override_sample_molecule_count']}",
                  file=sys.stderr)
            valid = False

        return valid

    def validate_and_set_optional_inputs(self, optional_inputs):
        """Helper method to validate and set optional inputs.

        Parameters
        ----------
        optional_inputs : dict
            Optional input files specified in the config file.

        Returns
        -------
        boolean
            True for valid optional inputs and False otherwise.

        """

        valid = True

        # Initialize list of all optional inputs
        self.optional_inputs['phased_vcf_file'] = None

        # If user hasn't specified any optional inputs, skip all validation steps.
        # This prevents errors when trying to access different elements of
        # optional_inputs.
        if optional_inputs is None:
            return valid

        # Check for presence of phased VCF file and validate existence.
        self.optional_inputs['phased_vcf_file'] = optional_inputs.get('phased_vcf_file', None)
        if self.optional_inputs['phased_vcf_file'] is not None:
            if not os.path.exists(self.optional_inputs['phased_vcf_file']) or \
               not os.path.isfile(self.optional_inputs['phased_vcf_file']):
                print(f"The phased vsf file specified under the optional inputs,"
                      f" {self.optional_inputs['phased_vcf_file']}, does not exist.",
                      file=sys.stderr)
                valid = False

        return valid

    def validate_and_set_sample_optional_inputs(self, sample_inputs):
        """Helper method to validate and set per sample optional inputs.

        Parameters
        ----------
        sample_inputs : dict
            Subject input parameters specified in the 'data' section of the
            config file.

        Returns
        -------
        boolean
            True for valid per sample optional inputs and False otherwise.

        """

        valid = True

        # TODO: incorporate BAM files into this optional input framework.

        # Still treat optional BAM files as a special case, handled by the
        # CAMPAREE controller class.
        available_optional_inputs = [
            # 'bam_file': None,
            'intron_quant',
            'gene_quant',
            'psi_quant',
            'allele_quant'
        ]

        for sample in self.samples:
            # Initialize dictionary  of optional files for current subject and
            # set all optional inputs to None.
            self.sample_optional_inputs[sample.sample_id] = {key: None for key in available_optional_inputs}
            # Skip to next sample if the current one has no optional inputs
            if sample_inputs[sample.sample_name]['optional_inputs'] is None:
                continue
            for optional_input, filename in sample_inputs[sample.sample_name]['optional_inputs'].items():
                # Treat optional BAM files as a special case, handled by the
                # CAMPAREE controller class.
                if optional_input == 'bam_file':
                    continue
                if optional_input in available_optional_inputs:
                    if filename is None:
                        # No filename given
                        print(f"No filename specified for optional input "
                              f"'{optional_input}', in config file for sample "
                              f"{sample.sample_name}.",
                              file=sys.stderr)
                        valid = False
                    elif not os.path.exists(filename) or not os.path.isfile(filename):
                        # Given file does not exist
                        print(f"File specified for optional input '{optional_input}'"
                              f" of sample {sample.sample_name},\n{filename}\n"
                              f"does not exist.",
                              file=sys.stderr)
                        valid = False
                    else:
                        self.sample_optional_inputs[sample.sample_id][optional_input] = filename
                else:
                    # Note, sort list of dictionary keys so this error message
                    # is consistent. The keys in default dictionaries are not
                    # maintained/returned in a sorted order.
                    print(f"Unrecognized optional input '{optional_input}', in "
                          f"config file for sample {sample.sample_name}.\n"
                          f"Must be one of the following: {','.join(sorted(available_optional_inputs))}.",
                          file=sys.stderr)
                    valid = False

        return valid


    def set_third_party_software(self):
        """
        Helper method to gather the names of all the 3rd party application files or directories and use them to set
        all the paths needed in the pipeline.  Since the third party software is shipped with this application, validation
        should not be necessary.  Software is identified generally by name and not specifically by filename since
        filenames may contain versioning and other artefacts.
        :return: the filenames for beagle, star, and kallisto, and the directory name for bowtie2.
        """

        beagle_filename = [filename for filename in os.listdir(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH)
                           if "beagle" in filename.lower()][0]
        self.beagle_file_path = os.path.join(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH, beagle_filename)

        star_filename = [filename for filename in os.listdir(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH)
                         if "STAR" in filename][0]
        self.star_file_path = os.path.join(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH, star_filename)

        kallisto_filename = [filename for filename in os.listdir(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH)
                              if "kallisto" in filename][0]
        self.kallisto_file_path = os.path.join(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH, kallisto_filename)

        bowtie2_dir_name = [filename for filename in os.listdir(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH)
                             if "bowtie2" in filename][0]
        self.bowtie2_dir_path = os.path.join(ExpressionPipeline.THIRD_PARTY_SOFTWARE_DIR_PATH, bowtie2_dir_name)

    def validate_and_set_resources(self, resources):
        """
        Since the resources are input file intensive, and since information about resource paths is found in the
        configuration file, this method validates that all needed resource information is complete, consistent and all
        input data is found.
        :param resources: dictionary containing resources from the configuration file
        :return: True if valid and False otherwise
        """
        # TODO a some point STAR will be in third party software and may require validation

        valid = True

        # Insure that all required resources keys are in place.  No point in continuing until this
        # problem is resolved.
        missing_resource_keys = [item
                                 for item in ExpressionPipeline.REQUIRED_RESOURCE_MAPPINGS
                                 if item not in resources]
        if missing_resource_keys:
            print(f"The following required mappings were not found under 'resources': "
                  f"{(',').join(missing_resource_keys)}", file=sys.stderr)
            return False

        # If user did not provide a path to the resources directory, use the
        # directory contained in the CAMPAREE install path.
        resources_directory_path = resources.get('directory_path', None)
        if not resources_directory_path:
            resources_directory_path = os.path.join(CAMPAREE_CONSTANTS.CAMPAREE_ROOT_DIR, "resources")
        elif not(os.path.exists(resources_directory_path) and os.path.isdir(resources_directory_path)):
            print(f"The given resources directory, {resources_directory_path}, must exist as a directory.",
                  file=sys.stderr)
            return False

        # Insure that the species model directory exists.  No point in continuing util this problem is
        # resolved.
        species_model_directory_path = os.path.join(resources_directory_path, resources['species_model'])
        if not(os.path.exists(species_model_directory_path) and os.path.isdir(species_model_directory_path)):
                print(f"The species model directory, {species_model_directory_path}, must exist as a directory",
                      file=sys.stderr)
                return False

        # Insure that the reference genome file path exists and points to a file.
        self.reference_genome_file_path = os.path.join(species_model_directory_path, resources['reference_genome_filename'])
        if not(os.path.exists(self.reference_genome_file_path) and os.path.isfile(self.reference_genome_file_path)):
            print(f"The reference genome file path, {self.reference_genome_file_path}, must exist as"
                  f" a file.", file=sys.stderr)
            valid = False

        # Insure that the chromosome ploidy file path exists and points to a file.
        self.chr_ploidy_file_path = os.path.join(species_model_directory_path, resources['chr_ploidy_filename'])
        if not(os.path.exists(self.chr_ploidy_file_path) and os.path.isfile(self.chr_ploidy_file_path)):
            print(f"The chr ploidy file path, {self.chr_ploidy_file_path} must exist as a file", file=sys.stderr)
            valid = False

        # Insure that the annotations file path exists and points to a file.
        self.annotation_file_path = os.path.join(species_model_directory_path, resources['annotation_filename'])
        if not(os.path.exists(self.annotation_file_path) and os.path.isfile(self.annotation_file_path)):
            print(f"The annotation file path, {self.annotation_file_path} must exist as a file", file=sys.stderr)
            valid = False

        # Insure that the star index directory exists as a directory.
        self.star_index_directory_path = \
            os.path.join(species_model_directory_path, resources['star_genome_index_directory_name'])
        if not (os.path.exists(self.star_index_directory_path) and os.path.isdir(self.star_index_directory_path)):
            print(f"The star index directory, {self.star_index_directory_path}, must exist as a directory",
                  file=sys.stderr)
            valid = False

        return valid

    def validate(self):
        valid = True
        reference_genome_chromosomes = self.reference_genome.keys()
        ploidy_chromosomes = set(self.chr_ploidy_data.keys())
        if not ploidy_chromosomes.issubset(reference_genome_chromosomes):
            missing_chromosomes = ' '.join(chrom for chrom in ploidy_chromosomes.difference(reference_genome_chromosomes))
            reference_chroms = ' '.join(chrom for chrom in reference_genome_chromosomes)
            print(f"The chromosome ploidy has chromosomes `{missing_chromosomes}` not found in the reference genome file", file=sys.stderr)
            print(f"The reference genome has chromosomes {reference_chroms}", file=sys.stderr)
            valid = False
        if not all([step.validate() for step in self.steps.values()]):
            valid = False
        return valid

    def execute(self):
        print("Execution of the Expression Pipeline Started...")

        seeds = self.generate_job_seeds()

        bam_files = {}
        for sample in self.samples:
            #Retrieve name of BAM file associated with this sample. This is either
            #the path to a user provided BAM file, or the default path the
            #GenomeAlignmentStep will use to store the alignment results for this
            #sample.
            genome_alignment = self.steps['GenomeAlignmentStep']
            bam_file = genome_alignment.get_genome_bam_path(sample)
            bam_files[sample.sample_id] = bam_file
            self.run_step(step_name='GenomeAlignmentStep',
                          sample=sample,
                          execute_args=[sample, self.star_index_directory_path,
                                        self.star_file_path],
                          cmd_line_args=[sample, self.star_index_directory_path,
                                         self.star_file_path])

        for sample in self.samples:
            bam_filename = bam_files[sample.sample_id]
            self.run_step(step_name='GenomeBamIndexStep',
                          sample=sample,
                          execute_args=[sample, bam_filename],
                          cmd_line_args=[sample, bam_filename],
                          dependency_list=[f"GenomeAlignmentStep.{sample.sample_id}"])

        for sample in self.samples:
            bam_filename = bam_files[sample.sample_id]
            seed = seeds[f"VariantsFinderStep.{sample.sample_id}"]
            self.run_step(step_name='VariantsFinderStep',
                          sample=sample,
                          execute_args=[sample, bam_filename, self.chr_ploidy_data,
                                        self.reference_genome, seed],
                          cmd_line_args=[sample, bam_filename, self.chr_ploidy_file_path,
                                         self.reference_genome_file_path, seed],
                          dependency_list=[f"GenomeBamIndexStep.{sample.sample_id}"])

        for sample in self.samples:
            if self.sample_optional_inputs[sample.sample_id]['intron_quant'] is None:
                output_directory = os.path.join(self.data_directory_path, f"sample{sample.sample_id}")
                bam_filename = bam_files[sample.sample_id]
                self.run_step(step_name='IntronQuantificationStep',
                              sample=sample,
                              execute_args=[bam_filename, output_directory, self.annotation_file_path],
                              cmd_line_args=[bam_filename, output_directory, self.annotation_file_path],
                              dependency_list=[f"GenomeBamIndexStep.{sample.sample_id}"])
            #TODO: do we need to depend upon the index being done? or just the alignment?
            #      I'm hypothesizing that some failures are being caused by indexing and quantification happening
            #      on the same BAM file at the same time, though I don't know why this would be a problem.

        seed = seeds["VariantsCompilationStep"]
        self.run_step(step_name='VariantsCompilationStep',
                      sample=None,
                      execute_args=[[sample.sample_id for sample in self.samples],
                                    self.chr_ploidy_data, self.reference_genome, seed],
                      cmd_line_args=[[sample.sample_id for sample in self.samples],
                                     self.chr_ploidy_file_path,
                                     self.reference_genome_file_path, seed],
                      dependency_list=[f"VariantsFinderStep.{sample.sample_id}" for sample in self.samples])

        phased_vcf_file = self.optional_inputs['phased_vcf_file']
        # If user did not provide phased vcf file
        if phased_vcf_file is None:
            seed = seeds["BeagleStep"]
            self.run_step(step_name='BeagleStep',
                          sample=None,
                          execute_args=[self.beagle_file_path, seed],
                          cmd_line_args=[self.beagle_file_path, seed],
                          dependency_list=[f"VariantsCompilationStep"])

        #TODO: We could load all of the steps in the entire pipeline into the queue
        #      and then just have the queue keep running until everything finishes.
        #      The only advantage of explicitly waiting here is that the user gets
        #      stdout indicating which stage is running for the pipeline.
        self.expression_pipeline_monitor.monitor_until_all_jobs_completed(queue_update_interval=10)

        for sample in self.samples:
            print(f"Processing sample{sample.sample_id} ({sample.sample_name})...")
            dep_list = None
            if phased_vcf_file is None:
                phased_vcf_file = os.path.join(self.data_directory_path,
                                                          CAMPAREE_CONSTANTS.BEAGLE_OUTPUT_FILENAME)
                dep_list = [f"BeagleStep"]
            self.run_step(step_name='GenomeBuilderStep',
                          sample=sample,
                          execute_args=[sample, phased_vcf_file, self.chr_ploidy_data,
                                        self.reference_genome],
                          cmd_line_args=[sample, phased_vcf_file, self.chr_ploidy_file_path,
                                         self.reference_genome_file_path],
                          dependency_list=dep_list)

            for suffix in [1, 2]:
                self.run_step(step_name='UpdateAnnotationForGenomeStep',
                              sample=sample,
                              execute_args=[sample, suffix, self.annotation_file_path,
                                            self.chr_ploidy_file_path],
                              cmd_line_args=[sample, suffix, self.annotation_file_path,
                                             self.chr_ploidy_file_path],
                              dependency_list=[f"GenomeBuilderStep.{sample.sample_id}"],
                              jobname_suffix=suffix)

                update_annot_path = os.path.join(self.data_directory_path, f"sample{sample.sample_id}",
                                                 CAMPAREE_CONSTANTS.UPDATEANNOT_OUTPUT_FILENAME_PATTERN.format(genome_name=suffix))
                parental_genome_path = os.path.join(self.data_directory_path, f"sample{sample.sample_id}",
                                                    CAMPAREE_CONSTANTS.GENOMEBUILDER_SEQUENCE_FILENAME_PATTERN.format(genome_name=suffix))
                self.run_step(step_name='TranscriptomeFastaPreparationStep',
                              sample=sample,
                              execute_args=[sample.sample_id, suffix, parental_genome_path,
                                            update_annot_path],
                              cmd_line_args=[sample.sample_id, suffix, parental_genome_path,
                                             update_annot_path],
                              dependency_list=[f"UpdateAnnotationForGenomeStep.{sample.sample_id}.{suffix}"],
                              jobname_suffix=suffix)

                tx_fasta_path = os.path.join(self.data_directory_path, f'sample{sample.sample_id}',
                                             CAMPAREE_CONSTANTS.TRANSCRIPTOME_FASTA_OUTPUT_FILENAME_PATTERN.format(genome_name=suffix))

                # Necessary for both gene and PSI quantification. Do not skip if
                # user provides optional input for only one of these distributions.
                if self.sample_optional_inputs[sample.sample_id]['gene_quant'] is None or \
                   self.sample_optional_inputs[sample.sample_id]['psi_quant'] is None:
                    self.run_step(step_name='KallistoIndexStep',
                                  sample=sample,
                                  execute_args=[sample.sample_id, suffix, self.kallisto_file_path,
                                                tx_fasta_path],
                                  cmd_line_args=[sample.sample_id, suffix, self.kallisto_file_path,
                                                 tx_fasta_path],
                                  dependency_list=[f"TranscriptomeFastaPreparationStep.{sample.sample_id}.{suffix}"],
                                  jobname_suffix=suffix)

                    self.run_step(step_name='KallistoQuantStep',
                                  sample=sample,
                                  execute_args=[sample, suffix, self.kallisto_file_path],
                                  cmd_line_args=[sample, suffix, self.kallisto_file_path],
                                  dependency_list=[f"KallistoIndexStep.{sample.sample_id}.{suffix}"],
                                  jobname_suffix=suffix)

                # Bowtie2 output only used for quantification of allelic imbalance.
                if self.sample_optional_inputs[sample.sample_id]['allele_quant'] is None:
                    self.run_step(step_name='Bowtie2IndexStep',
                                  sample=sample,
                                  execute_args=[sample.sample_id, suffix, self.bowtie2_dir_path,
                                                tx_fasta_path],
                                  cmd_line_args=[sample.sample_id, suffix, self.bowtie2_dir_path,
                                                 tx_fasta_path],
                                  dependency_list=[f"TranscriptomeFastaPreparationStep.{sample.sample_id}.{suffix}"],
                                  jobname_suffix=suffix)

                    self.run_step(step_name='Bowtie2AlignStep',
                                  sample=sample,
                                  execute_args=[sample, suffix, self.bowtie2_dir_path],
                                  cmd_line_args=[sample, suffix, self.bowtie2_dir_path],
                                  dependency_list=[f"Bowtie2IndexStep.{sample.sample_id}.{suffix}"],
                                  jobname_suffix=suffix)

            # Necessary for both gene and PSI quantification. Do not skip if
            # user provides optional input for only one of these distributions.
            if self.sample_optional_inputs[sample.sample_id]['gene_quant'] is None or \
               self.sample_optional_inputs[sample.sample_id]['psi_quant'] is None:
                # Use kallisto quantifications from the first parental genome to estimate
                # gene, transcript, and PSI quantifications for the simulated data.
                suffix=1
                kallisto_quant_path = os.path.join(self.data_directory_path, f'sample{sample.sample_id}',
                                                   CAMPAREE_CONSTANTS.KALLISTO_QUANT_DIR_PATTERN.format(genome_name=suffix),
                                                   CAMPAREE_CONSTANTS.KALLISTO_ABUNDANCE_FILENAME)
                update_annot_path = os.path.join(self.data_directory_path, f"sample{sample.sample_id}",
                                                 CAMPAREE_CONSTANTS.UPDATEANNOT_OUTPUT_FILENAME_PATTERN.format(genome_name=suffix))
                self.run_step(step_name='TranscriptGeneQuantificationStep',
                              sample=sample,
                              execute_args=[sample.sample_id, kallisto_quant_path, update_annot_path],
                              cmd_line_args=[sample.sample_id, kallisto_quant_path, update_annot_path],
                              dependency_list=[f"KallistoQuantStep.{sample.sample_id}.{suffix}"])

            if self.sample_optional_inputs[sample.sample_id]['allele_quant'] is None:
                genome_alignment_path = bam_files[sample.sample_id]
                update_annot_path_1 = os.path.join(self.data_directory_path, f"sample{sample.sample_id}",
                                                   CAMPAREE_CONSTANTS.UPDATEANNOT_OUTPUT_FILENAME_PATTERN.format(genome_name='1'))
                update_annot_path_2 = os.path.join(self.data_directory_path, f"sample{sample.sample_id}",
                                                   CAMPAREE_CONSTANTS.UPDATEANNOT_OUTPUT_FILENAME_PATTERN.format(genome_name='2'))
                tx_align_path_1 = os.path.join(self.data_directory_path, f'sample{sample.sample_id}',
                                               CAMPAREE_CONSTANTS.BOWTIE2_ALIGN_FILENAME_PATTERN.format(genome_name='1'))
                tx_align_path_2 = os.path.join(self.data_directory_path, f'sample{sample.sample_id}',
                                               CAMPAREE_CONSTANTS.BOWTIE2_ALIGN_FILENAME_PATTERN.format(genome_name='2'))
                self.run_step(step_name='AllelicImbalanceQuantificationStep',
                              sample=sample,
                              execute_args=[sample.sample_id, genome_alignment_path, update_annot_path_1,
                                            update_annot_path_2, tx_align_path_1, tx_align_path_2],
                              cmd_line_args=[sample.sample_id, genome_alignment_path, update_annot_path_1,
                                             update_annot_path_2, tx_align_path_1, tx_align_path_2],
                              dependency_list=[f"Bowtie2AlignStep.{sample.sample_id}.1",
                                               f"Bowtie2AlignStep.{sample.sample_id}.2"])

            seed = seeds[f"MoleculeMakerStep.{sample.sample_id}"]
            num_molecules_to_generate = sample.molecule_count
            # Use default or user-provided distribution files
            intron_quant_path = self.sample_optional_inputs[sample.sample_id]['intron_quant']
            gene_quant_path = self.sample_optional_inputs[sample.sample_id]['gene_quant']
            psi_quant_path = self.sample_optional_inputs[sample.sample_id]['psi_quant']
            allele_quant_path = self.sample_optional_inputs[sample.sample_id]['allele_quant']
            # These dependencies are only necessary if user provides optional
            # input for gene_quants, psi_quants, *and* allele_quants. This would
            # cause CAMPAREE to skip the steps used to generate those files. These
            # steps are dependent upon the TranscriptomeFastaPreparationStep and,
            # normally cause the pipeline to wait for the TranscriptomeFastaPreparationStep
            # to finish. If these steps that depend on TranscriptomeFastaPreparationStep
            # don't run, then the MoleculeMakerStep will start before the
            # TranscriptomeFastaPreparationStep has completed. This causes an
            # error because the molecule maker needs the transcriptome FASTA
            # produced by this step.
            dep_list = [f"TranscriptomeFastaPreparationStep.{sample.sample_id}.1",
                        f"TranscriptomeFastaPreparationStep.{sample.sample_id}.2"]
            if intron_quant_path is None:
                intron_quant_path = os.path.join(self.data_directory_path,
                                                 f'sample{sample.sample_id}',
                                                 CAMPAREE_CONSTANTS.INTRON_OUTPUT_FILENAME)
            if gene_quant_path is None:
                gene_quant_path = os.path.join(self.data_directory_path,
                                                 f'sample{sample.sample_id}',
                                                 CAMPAREE_CONSTANTS.TXQUANT_OUTPUT_GENE_FILENAME)
                dep_list.append(f"TranscriptGeneQuantificationStep.{sample.sample_id}")
            if psi_quant_path is None:
                psi_quant_path = os.path.join(self.data_directory_path,
                                                 f'sample{sample.sample_id}',
                                                 CAMPAREE_CONSTANTS.TXQUANT_OUTPUT_PSI_FILENAME)
                # Dependency already added if user did not provide gene quant.
                if f"TranscriptGeneQuantificationStep.{sample.sample_id}" not in dep_list:
                    dep_list.append(f"TranscriptGeneQuantificationStep.{sample.sample_id}")
            if allele_quant_path is None:
                allele_quant_path = os.path.join(self.data_directory_path,
                                                 f'sample{sample.sample_id}',
                                                 CAMPAREE_CONSTANTS.ALLELIC_IMBALANCE_OUTPUT_FILENAME)
                dep_list.append(f"AllelicImbalanceQuantificationStep.{sample.sample_id}")

            # If no molecule count specified for this sample, use the default count.
            if not num_molecules_to_generate or self.override_sample_molecule_count:
                num_molecules_to_generate = self.default_molecule_count
            self.run_step(step_name='MoleculeMakerStep',
                          sample=sample,
                          execute_args=[sample, intron_quant_path, gene_quant_path,
                                        psi_quant_path, allele_quant_path,
                                        self.output_type, num_molecules_to_generate, seed],
                          cmd_line_args=[sample, intron_quant_path, gene_quant_path,
                                         psi_quant_path, allele_quant_path,
                                         self.output_type, num_molecules_to_generate, seed],
                          dependency_list=dep_list)

        self.expression_pipeline_monitor.monitor_until_all_jobs_completed(queue_update_interval=10)

        print("Execution of the Expression Pipeline Ended")

    def generate_job_seeds(self):
        """
        Generate one seed per job that needs a seed, returns a dictionary mapping
        job names to seeds

        We generate seeds for each job since they run on separate nodes of the cluster, potentially
        and so do not simply share Numpy seeds. We generate them all ahead of time so that if jobs need
        to be restart, they can reuse the same seed.
        """
        seeds = {}
        # Seeds for jobs that don't run per sample
        for job in ["VariantsCompilationStep", "BeagleStep"]:
            seeds[job] = numpy.random.randint(MAX_SEED)
        # Seeds for jobs that are run per sample
        for job in ["VariantsFinderStep", "MoleculeMakerStep"]:
            for sample in self.samples:
                seeds[f"{job}.{sample.sample_id}"] = numpy.random.randint(MAX_SEED)
        return seeds

    def run_step(self, step_name, sample, execute_args, cmd_line_args, dependency_list=None,
                 jobname_suffix=None):
        """
        Helper function that runs the given step, with the given parameters. It
        wraps submission of the step to the scheduler/job monitor.

        Parameters
        ----------
        step_name : string
            Name of the CAMPAREE step to run. It should be in the list of steps
            stored in the steps dictionary.
        sample : Sample
            Sample to run through the step. For steps that aren't associated with
            specific samples, set this to None.
        execute_args : list
            List of positional paramteres to pass to the execute() method for the
            given step.
        cmd_line_args : list
            List of positional parameters to pass to the get_commandline_call()
            method for the given step.
        dependency_list : list
            List of job names (if any) the current step depends on. Default: None.
        jobname_suffix : string
            Suffix to add to job submission ID. Default: None.

        """
        if step_name not in list(self.steps.keys()):
            raise CampareeException(f"{step_name} not in the list of loaded steps (see config file).")

        step_class = self.steps[step_name]

        # Check if the current step has an overrides for scheduler parameters
        scheduler_num_processors = None
        scheduler_memory_in_mb = None
        if self.__step_scheduler_param_overrides[step_name]:
            scheduler_num_processors = self.__step_scheduler_param_overrides[step_name].get('num_processors',None)
            scheduler_memory_in_mb = self.__step_scheduler_param_overrides[step_name].get('memory_in_mb',None)

        # Check if there are any scheduler submission arguments to add
        scheduler_submission_args = ""
        if self.scheduler_default_params['default_submission_args']:
            scheduler_submission_args = self.scheduler_default_params['default_submission_args']

        stdout_log = os.path.join(step_class.log_directory_path,
                                  f"sample{sample.sample_id}" if sample else "",
                                  f"{step_name}{f'.{jobname_suffix}' if jobname_suffix else ''}.{self.scheduler_mode}.out")
        stderr_log = os.path.join(step_class.log_directory_path,
                                  f"sample{sample.sample_id}" if sample else "",
                                  f"{step_name}{f'.{jobname_suffix}' if jobname_suffix else ''}.{self.scheduler_mode}.err")
        scheduler_job_name = (f"{step_name}{f'.sample{sample.sample_id}_{sample.sample_name}' if sample else ''}"
                              f"{f'.{jobname_suffix}' if jobname_suffix else ''}")
        scheduler_args = {'job_name' : scheduler_job_name,
                          'stdout_logfile' : stdout_log,
                          'stderr_logfile' : stderr_log,
                          'num_processors' : scheduler_num_processors,
                          'memory_in_mb' : scheduler_memory_in_mb,
                          'additional_args' : scheduler_submission_args
                         }
        command = step_class.get_commandline_call(*cmd_line_args)
        validation_attributes = step_class.get_validation_attributes(*cmd_line_args)
        output_directory = os.path.join(step_class.data_directory_path,
                                        f"sample{sample.sample_id}" if sample else "")
        self.expression_pipeline_monitor.submit_new_job(job_id=f"{step_name}{f'.{sample.sample_id}' if sample else ''}"
                                                               f"{f'.{jobname_suffix}' if jobname_suffix else ''}",
                                                        job_command=command,
                                                        sample=sample,
                                                        step_name=step_name,
                                                        scheduler_arguments=scheduler_args,
                                                        validation_attributes=validation_attributes,
                                                        output_directory_path=output_directory,
                                                        system_id=None,
                                                        dependency_list=dependency_list)

    @staticmethod
    def main(configuration, scheduler_mode, output_directory_path, input_samples):
        pipeline = ExpressionPipeline(configuration, scheduler_mode, output_directory_path, input_samples)
        if not pipeline.validate():
            raise CampareeValidationException("Expression Pipeline Validation Failed.  "
                                              "Consult the standard error file for details.")
        pipeline.execute()


class CampareeValidationException(CampareeException):
    pass
