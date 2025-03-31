# This script creates a baby.config.yaml file in the local path
# which contains an example configuration for CAMPAREE and is used
# to test the installation for sucess.

import pathlib
from camparee.camparee_constants import CAMPAREE_CONSTANTS

def main():
    config_dir = pathlib.Path(CAMPAREE_CONSTANTS.CAMPAREE_ROOT_DIR) / "config"
    baby_config = config_dir / "baby.config.yaml"
    new_file = pathlib.Path("baby.config.yaml")

    new_file.write_text(baby_config.read_text())
    print(f"Created {str(new_file)}")
