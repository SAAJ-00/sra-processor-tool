# setup.py
from setuptools import setup, find_packages

setup(
    name="sra_processor",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        'pathlib',
        'argparse',
    ],
    python_requires='>=3.6',
    entry_points={
        'console_scripts': [
            'sra-processor=sra_processor.cli:main',
        ],
    },
)