from setuptools import setup, find_packages

setup(
    name="blasch",
    version="0.1.0",
    description="Chimera detection and recovery using BLAST",
    author="Ali Hakimzadeh",
    packages=find_packages(),
    install_requires=[
        "biopython",
    ],
    entry_points={
        'console_scripts': [
            'blasch=blasch.main:main',
        ],
    },
    python_requires=">=3.6",
)