Running the Pipeline
********************

This section explores various ways to actually run the pipeline. Although some of these options are described in Command Line options, this will delve into more details

Luigi
-----

As mentioned before, the dependency management of the pipeline is handled through a python package called luigi. However, currently the only options available when running luigi are adjusting the number of workers. You can read more about workers `here <https://luigi.readthedocs.io/en/stable/api/luigi.worker.html>`_

Softwares
---------

There are two options for specifying where the softwares' folder location will be. 

Using Slurm Jobs
----------------

Tmux
----