#!/bin/python
from setuptools import find_packages, setup

setup(
    name="PID_pendulum",
    version="0.0.1",
    author="Philipp Schuette",
    author_email="p.schuette@online.e",
    description="A PID pendulum controller.",
    url="https://github.com/PhilippSchuette/PID_pendulum",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GPL-3.0 License",
        "Operating System :: OS Independent",
    ]
)
