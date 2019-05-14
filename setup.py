#!/bin/python
from setuptools import find_packages, setup

setup(
    name="PID_pendulum",
    version="0.0.3",
    author="Philipp Schuette",
    author_email="p.schuette@online.de",
    description="A PID pendulum controller.",
    url="https://github.com/PhilippSchuette/PID_pendulum",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        "matplotlib",
        "numpy"
    ]
)
