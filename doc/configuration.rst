Command Line Options
********************

This section explores the different options that are possible when running the tool. Before running the tool, you can first run the command ``mitopipeline -d`` in order to download all of the dependencies necessary for the steps in your pipeline if you do not specify a tools directory. 

Downloading Software
--------------------

.. csv-table::
    :header: "Option", "Description"
    :widths: 20, 30

    "``-d`` or ``--download``", "| **REQUIRED**. 
    | Option to download softwares for steps not listed in -r."
    "``-r`` or ``--remove``", "| Specifies steps that won't be run in the pipeline. 
    | **Input**: <name of steps> 
    | **Default**: None" 

Running Pipeline
------------------

.. csv-table::
    :header: "Option", "Description"
    :widths: 20, 30

    "``-s`` or ``--start``", "| **REQUIRED**. 
    | Specifies the directory that contains the files to be run on. 
    | **Input**: <path/to/startdirectory>"
    "``-g``, ``--genomes``", "| **REQUIRED**.
    | Specifies location of the folder that contains the human reference genomes. 
    | The required genomes for this pipeline are: 
    | - hg38-nochr.fa (the human genome without the mitochondrial genome)
    | - hg38.fa (the entire human genome)
    | - rCRS.fa (the human mitochondrial genome). 
    | The file names must match the ones listed here so that the steps are able to find the files. You can also download these human references genomes from the `UCSC genome browser <http://hgdownload.cse.ucsc.edu/downloads.html#human>`_.
    | **Input:** <path/to/REFs/directory>"
    "``-t``, ``---tools``", "| **REQUIRED**.
    | Specifies the location of the folder that contains all of the 3rd party software. The executables must be in a folder with the same name as the software, i.e. gatk should have its default executable within </path/to/tools/gatk> .
    | **Input**: <path/to/tools/directory>
    | **Default**: <path/to/mitopipeline/tools/directory>"
    "``-o`` or ``--output``", "| Specifies where to store the output of the pipeline results. Will create folders for each step within the output directory.
    | **Input**: <path/to/output/directory>
    | **Default**: current directory"
    "``-l`` or ``--slurm``", "| Use the slurm workload manager to submit batch jobs for each step of the pipeline
    | **Default**: False"
    "``-w`` or ``--workers``", "| Specifies the number of workers to use when running the pipeline.
    | **Input**: integer
    | **Default**: 1"
    "``-r`` or ``--remove``", "| Specifies steps that won't be run in the pipeline. The names of the steps are:
    | ['extractmito', 'splitgap', 'clipping', 'removenumts', 'gatk', 'snpeff', 'annovar']
    | **Input**: <name of steps> 
    | **Default**: None"