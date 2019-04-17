Getting Started
---------------
This project can be downloaded through ``pip install mitopipeline`` to install the latest version of mitopipeline. `Documentation <https://mitopipeline.readthedocs.io/en/latest/>`_ is hosted through readthedocs.

Users will have to download samtools and bwa, ``make`` them, and include the path to the executable in their $PATH variable or copy them to their ``bin`` folder so that they can be run from the command line. Download links for `samtools  <http://www.htslib.org/download/>`_ and `bwa <https://sourceforge.net/projects/bio-bwa/>`_. In addition, 3rd party software packages that they wish to use on their own (such as GATK, ANNOVAR, etc.) need to be downloaded before being able to use them in the pipeline.

Purpose
-------
The purpose of this tool is to provide a useful pipeline for genetic researchers working with the mitochondrial genome. 

In genetics research, raw genome data requires many different softwares to be run first to extract useful information such as variant information before it can be analyzed. Thus, we want to provide a pipeline that allows users to streamline their processing automatically for them. 

The steps included within this pipeline include extracting mitochondrial genome from human genome, clipping, splitting the gap for Complete Genomics data, alignment into NuMT removals, and other software packages such as GATK, SNPEFF, ANNOVAR.

Goals
-----
This tool should handle all dependency management in-between each step for the user from start to finish with the end result being processed data that the user can begin analyzing immediately. Specifically, this tool should:

- Handle all dependencies within the pipeline from start to finish

- Provide built in steps such as alignment, clipping, removal of NuMTs, and splitting the gap for Complete Genomics vcf files

- Use 3rd party software packages like GATK, ANNOVAR, SNPEFF, etc.

- Accurately log successes and failures that occurred while running the pipeline

Credits
----------------
The dependency management is controlled through a python package called LuigiMany of the steps in the pipeline were adapted from scripts written by Sneha Grandhi. Check out her `Github  <https://github.com/sneha-grandhi/>`_, `LinkedIn <https://www.linkedin.com/in/sneha-grandhi-phd-0165aa58/>`_, or contact her through her  email: sneha_grandhi@hms.harvard.edu.

Table of Contents
-----------------