#!/bin/python
from setuptools import find_packages, setup

version = "0.0.10"
description = "Visit https://github.com/PhilippSchuette/PID_pendulum for" \
    " additional information."

setup(
    name="PID_pendulum",
    version=version,
    author="Philipp Schuette",
    author_email="p.schuette@online.de",
    description="A PID pendulum controller.",
    long_description=description,
    long_description_content_type="text/plain",
    url="https://github.com/PhilippSchuette/PID_pendulum",
    download_url="https://pypi.python.org/pypi/PID-pendulum",
    license="GPLv3",
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
