Running the Pipeline
********************

This section explores details about running the pipeline such as the tools it uses, its requirements, and some suggestions on how to run it. Although some of these options are described in Command Line options, this will delve into more details.

Luigi
-----

As mentioned before, the dependency management of the pipeline is handled through a python package called `luigi <https://github.com/spotify/luigi>`_. However, currently the only options available when running luigi are adjusting the number of workers. You can read more about workers `here <https://luigi.readthedocs.io/en/stable/api/luigi.worker.html>`_

Luigi handles dependency management through class objects called Targets and Tasks. Targets are a Task’s output. A Task is run after a required Task is complete and also outputs a Target for a next Task to be run. For example, the workflow for two tasks running on a database can be shown like this:	


.. figure:: https://raw.githubusercontent.com/timmykuo/mitopipeline/master/doc/luigi_tasks_targets.png


In this diagram, the first task takes in the data from the database as input and outputs a target. The target is then input to the next task to be run. In order for the second task to be run, it “requires” the first task to be finished first. This is tracked through the existence of the first task’s output (the target). Once it sees the target in the output, the 2nd task will begin  running. The advantage of such a design is its asynchronous processes. Since the time for each individual file may be different for the same task, having a worker that looks solely for the output target allows for the multiple tasks to be run at the same time.

Softwares
---------

As descrbied in the pipeine steps section, all of the steps have some software requirements in order to be run. There are two options for getting the softwares necessary. 

The first choice is to use the command line option ``-d``. For example, the command

.. code:: console

    $ mitopipeline -d -r annovar snpeff

will download all the necessary software into mitopipeline's tool's directory for all steps except for annovar and snpeff. You can then use the mitopipeline normally.

The second choice is to specify a directory that has all the necessary softwares downloaded. Keep in mind that mitopipeline will check for the naming convention of the software's folder that contains its executable as the same name as the step i.e. 'gatk' step will look for a folder called 'gatk' within the specified directory for a gatk.jar executable. 

The pipeline will also check for if an executable is available from the command line. For example, if ``samtools`` and ``bwa`` can be executed on the command line (in your $PATH) and is not in the tools directory, that is acceptable.

Using Slurm Jobs
----------------

Some servers have the `slurm workload manager <https://slurm.schedmd.com/overview.html>`_ set up on their system. If you are using such a server, an available option is to use the option ``-l``. This will submit slurm jobs for each step of the pipeline for each file and save the files in a folder within the specified -out directory.

Tmux
----

Currently, luigi's scheduler is not implemented within this tool and only uses its local scheduler (read in luigi's docs). Thus, it requires that whatever process that is running mitopipeline to be continually running. One way to do this is to run it on a server using a tmux session. You can read more about tmux here.

Once tmux is downloaded, you can start a new tmux session by typing ``tmux`` into your command line. Then, after beginning the pipeline through the ``mitopipeline`` command, you can exit the session by pressing ``ctrl+b`` and then ``d``. This will detach the current tmux session from your terminal.

In order to reenter your tmux session, you can type in ``tmux ls`` in order to list all of your sessions and then ``tmux a -t <your-session-id>`` to re-enter that tmux session where your mitopipeline is running.