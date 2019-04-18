Command Line Options
********************

This section explores the different options that are possible when running the tool. Before running the tool, you can first run the command ``mitopipeline -d`` in order to download all of the dependencies necessary for the steps in your pipeline if you do not specify a tools directory. 

Downloading Software
--------------------

.. csv-table::
    :header: "Option", "Description"

    "``-d`` or ``--download``", "| **REQUIRED**. 
    | Option to download softwares for steps not listed in -r."
    "``-r`` or ``--remove``", "| Specifies steps that won't be run in the pipeline. 
    | **Input**: <name of steps> 
    | **Default**: None" 

Running Pipeline
------------------

.. csv-table::
    :header: "Option", "Description"

    "``-s`` or ``--start``", "| **REQUIRED**. 
    | Specifies the directory that contains the files to be run on. 
    | **Input**: <path/to/startdirectory>"
    "``-t``, ``---tools``", "| Specifies the location of the folder that contains all of the 3rd party software.
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
    "``-r`` or ``--remove``", "| Specifies steps that won't be run in the pipeline. 
    | **Input**: <name of steps> 
    | **Default**: None"