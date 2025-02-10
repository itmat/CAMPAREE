#!/bin/bash

wd="$(pwd)"

# Handle options
options=$(getopt -o i: --long id: -- "$@")
eval set -- "${options}"
while true; do
    case "${1}" in
        -i|--id)
            id="${2}"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Invalid option: ${1}"
            exit 1
            ;;
    esac
done

camparee_resources_path=$(find "${SP_JOB_RESOURCES_DIR}/camparee/" -name "GRCh38_gencode_v46_primary_assembly.tgz")

mkdir resources

tar -xzf "${camparee_resources_path}" -C resources

# Finding trimmed FASTQ files
echo "Finding trimmed FASTQ files"
fastq_dir="$(dirname $(find . -type f -regextype posix-extended -iregex '.*_1\.trimmed.*\.fastq(\.gz)?$'))"
fastq_dir="${fastq_dir#*/}"
trimmed_file_1="$(basename $(find . -type f -regextype posix-extended -iregex '.*_1\.trimmed.*\.fastq(\.gz)?$'))"
trimmed_file_2="$(basename $(find . -type f -regextype posix-extended -iregex '.*_2\.trimmed.*\.fastq(\.gz)?$'))"
#VALIDATE_FILE "${trimmed_file_1}"

python3 create_config.py GRCh38_gencode_v46_primary_assembly.config.yaml \
"${wd}" "${id}" "${fastq_dir}" "${trimmed_file_1}" "${trimmed_file_2}" "${id}.config.yaml"

mkdir results
config_file="${wd}/${id}.config.yaml"
cd /CAMPAREE && poetry run bin/run_camparee.py -c "${config_file}"

cd "${wd}"
mv results/run_1/CAMPAREE/data/sample1/molecule_file.txt ${id}_molecule_file.txt
