import os

from setuptools import setup

setup(
    name="MitoPipeline",
    version="0.1",
    packages=find_packages(),
    scripts=[],

    #add samtools, bwa, etc. here
    install_requires=[],

    package_data={}

    #metadata
    author="Timothy Kuo",
    author_email="timykuo@gmail.com",
    description="Pipeline Creator",
    license="MIT",
    keywords="genetics research data pipeline dependency management",
    url="",
    project_urls={
        "Bug Tracker": "",
        "Documentation": "",
        "Source Code": "",
    }
)