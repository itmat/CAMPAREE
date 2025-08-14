# This script creates a baby.config.yaml and all the files it points to in the local path
# which contains an example configuration for CAMPAREE and is used
# to test the installation for sucess.

import pathlib
import shutil
from camparee.camparee_constants import CAMPAREE_CONSTANTS

def main():
    config_dir = (pathlib.Path(CAMPAREE_CONSTANTS.CAMPAREE_ROOT_DIR) / "config").resolve()
    baby_config = config_dir / "baby.config.yaml"
    test_data_dir = (pathlib.Path(CAMPAREE_CONSTANTS.CAMPAREE_ROOT_DIR) / "test_data").resolve()
    resources = (pathlib.Path(CAMPAREE_CONSTANTS.CAMPAREE_ROOT_DIR) / "resources").resolve()

    new_config = pathlib.Path("baby.config.yaml")

    new_config.write_text(baby_config.read_text())
    print(f"Created {str(new_config)}")

    shutil.copytree(
            test_data_dir,
            "test_data",
            dirs_exist_ok=True
        )
    print("Created test_data/")
    shutil.copytree(
            resources,
            "resources",
            dirs_exist_ok=True
        )
    print("Created resources/")
