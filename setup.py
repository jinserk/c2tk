import setuptools
import distutils.text_file as txt
from pathlib import Path
from typing import List

def _parse_requirements(filename: str) -> List[str]:
    """Return requirements from requirements file."""
    # Ref: https://stackoverflow.com/a/42033122/
    return txt.TextFile(filename=str(Path(__file__).with_name(filename))).readlines()

setuptools.setup(
    install_requires=_parse_requirements('requirements.txt')
)
