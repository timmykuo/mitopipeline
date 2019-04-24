Getting Started
---------------
This project can be downloaded through ``pip install mitopipeline`` to install the latest version of mitopipeline. Or, you can directly clone the repository through ``git clone https://github.com/timmykuo/mitopipeline.git``, navigate to the cloned repository, and run ``python3 setup.py install``. `Documentation <https://mitopipeline.readthedocs.io/en/latest/>`_ is hosted through readthedocs.

Users will have to download samtools and bwa, ``make`` them, and include the path to the executable in their $PATH variable or copy them to their ``bin`` folder so that they can be run from the command line. Download links for `samtools  <http://www.htslib.org/download/>`_ and `bwa <https://sourceforge.net/projects/bio-bwa/>`_. In addition, 3rd party software packages that they wish to use on their own (such as GATK, ANNOVAR, etc.) need to be downloaded before being able to use them in the pipeline.

Purpose
-------
In genetics research, researchers often require raw genome data to be run on different softwares first to extract useful information such as variant information before it can be analyzed. Thus, we want to provide a pipeline that allows users to streamline their processing automatically for them. 

The steps included within this pipeline include extracting mitochondrial genome from human genome, clipping, splitting the gap for Complete Genomics data, alignment into NuMT removals, and other software packages such as GATK, SNPEFF, ANNOVAR.

Pipeline
--------
.. figure:: https://raw.githubusercontent.com/timmykuo/mitopipeline/master/doc/pipeline.png

The steps shown in this figure are described further in the Pipeline Steps page. Each step can be omitted if it's unncessary for your pipeline and can also be replaced with your own scripts.

Credits
----------------
The dependency management is controlled through a python package called `Luigi <https://github.com/spotify/luigi/>`_. In addition, many of the steps in the pipeline were adapted from scripts written by Sneha Grandhi. Check out her `Github  <https://github.com/sneha-grandhi/>`_, `LinkedIn <https://www.linkedin.com/in/sneha-grandhi-phd-0165aa58/>`_, or contact her through her  email: sneha_grandhi@hms.harvard.edu.

Table of Contents
-----------------