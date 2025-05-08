"""
This module defines how the data directory for the library is determined.

By default, it uses a subdirectory named `data` located relative to the current file.
However, this can now be overridden by user configuration.
"""

import os
import platform
from pathlib import Path
import yaml


def get_downscaled_data_dir(override: str = None) -> str:
    """
    Determines the directory path for data storage based on various sources of configuration.

    This function checks for the data directory in the following order of precedence:

    1. If an override path is provided as an argument, it is returned.
    2. If a config file exists in standard locations, the `data_dir` from the file is returned.

       - On Windows, it checks the `APPDATA` directory.
       - On other platforms, it checks `~/.config/tcr/config.yaml`.

    3. If none of the above sources provide a path, a default `DATA_DIR` constant is returned.

    Args:
        override (str, optional): A user-specified directory path to override all other sources.
        Defaults to None.

    Returns:
        str: The resolved data directory path.
    """
    if override:
        return override

    config_paths = []

    # Cross-platform config locations
    if platform.system() == "Windows":
        appdata = os.getenv("APPDATA")  # Typically: C:\Users\<User>\AppData\Roaming
        if appdata:
            config_paths.append(Path(appdata) / "tcr" / "config.yaml")
    else:
        config_paths.append(Path.home() / ".config" / "tcr" / "config.yaml")

    for config_path in config_paths:
        if config_path.is_file():
            with open(config_path, "r", encoding="utf-8") as f:
                config = yaml.safe_load(f)
                if "data_dir" in config:
                    return config["data_dir"]

    return BASE_DATA_DIR


# Default data directory relative to this file
BASE_DATA_DIR = os.path.join(os.path.split(__file__)[0], 'data')

# Downscaled data directory from yaml config
# This is a separate directory for downscaled data, which may be different from the main directory
DOWNSCALED_DATA_DIR = get_downscaled_data_dir()
