# This script uses the STAR executable packaged with CAMPAREE to create
# a STAR index for the baby genome. This index is created in the
# resources/baby_genome.mm10/star_index.genome/ directory in the CAMPAREE install
# path.

import sys
import pathlib
import subprocess
from camparee.camparee_constants import CAMPAREE_CONSTANTS

def main():
    resources_dir = pathlib.Path("resources")
    if not resources_dir.exists() and resources_dir.is_dir():
        print("Expected 'resources/' directory to be present. Did you run `create_camparee_test_files` command firt?")
        sys.exit()
    star_path = pathlib.Path(CAMPAREE_CONSTANTS.CAMPAREE_ROOT_DIR) / "third_party_software" / "STAR"

    # Create directory for babY_gneome STAR index (if it doesn't already exist)
    star_index_dir = resources_dir / "baby_genome.mm10/star_index.genome"
    star_index_dir.mkdir(exist_ok=True)


    #Run STAR to create the index (if it's destination directory was created correctly)
    subprocess.run(
        [
            star_path, 
            "--runMode", "genomeGenerate",
            "--runThreadN", "1",
            "--genomeDir", str(star_index_dir),
            "--outFileNamePrefix", f"{str(star_index_dir)}/",
            "--genomeFastaFiles", str(resources_dir / "baby_genome.mm10" / "baby_genome.mm10.oneline_seqs.fa"),
        ],
        check = True
    )
