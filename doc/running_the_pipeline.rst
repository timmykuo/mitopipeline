Running the Pipeline
********************

This section explores details about running the pipeline such as the tools it uses, its requirements, and some suggestions on how to run it. Although some of these options are described in Command Line options, this will delve into more details.

Required Arguments
------------------

| **Start Directory:** 
| The start directory is necessary in order to specify where all the files that need to be run are. One thing to note is that the pipeline use's the '.' from the file's extension to parse the file name. Thus, all files within the start directory can only have one period within its name that specifies the extension, i.e. FILENAME.bam. In addition, .bai index files are allowed to be in the same directory as the .bam files.

| **Tools Directory:** 
| In order to run 3rd party softwares like seqtk, gatk, snpeff, and annovar, mitopipeline requires you to specify where these tools are kept. Each software package should have a folder within the tools directory that contains its executables. For example, gatk's executable should be path/to/tools/gatk/GenomeAnalysisTK.jar, snpeff's should be path/to/tools/snpeff/snpEff.jar and annovar's should be path/to/tools/annovar/annovar-executables. One choice is to first use mitopipeline's download function if you don't want to download these yourself.

| **Reference Genomes:**
| With steps such as GATK and RemoveNuMTs, it's necessary to have human reference genomes to align to. The required genomes for this pipeline are `hg38-nochr.fa <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/>`_ (the GRCh38/hg38 version human genome without the mitochondrial genome), `hg38.fa <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/>`_ (the GRCh38/hg38 version of the human mitochondrial genome), and `rCRS-MT.fa <http://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz>`_ (the GRCh38/hg38 version of the human mitochondrial genome). The file names must be changed to match the ones listed here so that the steps are able to find the files. You can read more about the human references genomes from the `UCSC genome browser <http://hgdownload.cse.ucsc.edu/downloads.html#human>`_. Since the files are too large to be downloaded along with the mitopipeline package, the user must specify the path to the reference genomes after downloading.

Python Modules on Case's HPC
-----------------------------

On an online server, users typically do not have permission to download softwares and modules directly due to lack of administrative privileges. In this case, users will need to set up python modules that contains the location of mitopipeline. Here's an example of how to do this on Case Western's HPC server that uses the slurm workload mananger and lmod. For more details, you can refer to the Case Western HPC website for more detailed `instructions <https://sites.google.com/a/case.edu/hpcc/hpc-tutorials/installing-local-python-modules>`_ on how to use pip to install and use python modules:

First, to install the module, we need to first set up the enviromental module $PYTHONUSERBASE.

.. code:: bash
    
    export PYTHONUSERBASE=$HOME/.usr/local/python/3.5.1
    pip install --user mitopipeline

This will install mitopipline into the $PYTHONUSERBASE directory. To use the installed module (since this tool contains binaries), we need to include it in the path. For our local directory of packages, we will create a module file that will set those variables for us. We will call the module "python-modules" and we will set the version of the module to the version of Python we are using. First, we will create the directory (just once).

.. code:: bash

    PYTHONMODULES=$HOME/.usr/local/share/modulefiles/python-modules
    mkdir -p $PYTHONMODULES
    cd $PYTHONMODULES

Then, we will create a file called 3.5.1-gcc.lua in this directory that contains the following contet:

.. code:: bash

    -- This is Lua module file for our local Python
    -- modules.
    -- To use it just run
    --    module load python-modules/3.5.1-gcc
    --
    --

    load("intel/17", "openmpi/2.0.1","python/3.5.1")

    pushenv("PYTHONUSERBASE",pathJoin(os.getenv("HOME"),".usr/local/python/3.5.1"))
    prepend_path("PATH",pathJoin(os.getenv("HOME"),".usr/local/python/3.5.1/bin))

To use this module, just run

..code:: bash
    module load python-modules/3.5.1-gcc

From here, if you used mitopipeline's "-d" command line option to download tools for you, the tools directory can now be referenced through: " /home/<case_ID>/.usr/local/python/3.5.1/lib/python3.5/site-packages/mitopipeline/tools/". Note: If you are using ANNOVAR, you must move the entire folder over to "~/.local/lib/python3.5/site-packages/mitopipeline" after downloading yourself, because ANNOVAR requires user registration before downloading (so it's unavailable through "-d"). Of course, you can also specify your own tools directory.

Luigi
-----

As mentioned before, the dependency management of the pipeline is handled through a python package called `luigi <https://github.com/spotify/luigi>`_. However, currently the only options available when running luigi are adjusting the number of workers. You can read more about workers `here <https://luigi.readthedocs.io/en/stable/api/luigi.worker.html>`_

Luigi handles dependency management through class objects called Targets and Tasks. Targets are a Task’s output. A Task is run after a required Task is complete and also outputs a Target for a next Task to be run. For example, the workflow for two tasks running on a database can be shown like this:	


.. figure:: https://raw.githubusercontent.com/timmykuo/mitopipeline/master/doc/luigi_tasks_targets.png


In this diagram, the first task takes in the data from the database as input and outputs a target. The target is then input to the next task to be run. In order for the second task to be run, it “requires” the first task to be finished first. This is tracked through the existence of the first task’s output (the target). Once it sees the target in the output, the 2nd task will begin  running. The advantage of such a design is its asynchronous processes. Since the time for each individual file may be different for the same task, having a worker that looks solely for the output target allows for the multiple tasks to be run at the same time.

Softwares
---------

As described in the pipeine steps section, all of the steps have some software requirements in order to be run. There are two options for getting the softwares necessary. 

The first choice is to use the command line option ``-d``. For example, the command

.. code:: console

    $ mitopipeline -d -r annovar snpeff

will download all the necessary software into mitopipeline's tool's directory for all steps except for annovar and snpeff. You can then use the mitopipeline normally without specifying the tools directory.

The second choice is to specify a directory that has all the necessary softwares downloaded. This is only necessary only for the step softwares, including seqtk, GATK, SNPEFF, and ANNOVAR. Keep in mind that mitopipeline will check for the naming convention of the software's folder that contains its executable as the same name as the step i.e. 'gatk' step will look for a folder called 'gatk' within the specified directory for its executable. 

A number of softwares are necessary to be run on the command line as they are called directly through the bash scripts. In particular, 'samtools' and 'bwa' need to be able to be executed through the comand line. On MacOSX/Linux, this can be achieved by either copying the executable to your ``/usr/local/bin`` folder or adding the folder of your executable to your $PATH variable. You can read more about each step's required softwares on the Pipeline Steps page.

Using Slurm Jobs
----------------

Some servers have the `slurm workload manager <https://slurm.schedmd.com/overview.html>`_ set up on their system. If you are using such a server, an available option is to use the option ``-l``. This will submit slurm jobs for each step of the pipeline for each file and save the files in a folder within the specified -out directory.

Tmux
----

Currently, luigi's scheduler is not implemented within this tool and only uses its local scheduler (read in luigi's docs). Thus, it requires that whatever process that is running mitopipeline to be continually running. One way to do this is to run it on a server using a tmux session. You can read more about tmux here.

Once tmux is downloaded, you can start a new tmux session by typing ``tmux`` into your command line. Then, after beginning the pipeline through the ``mitopipeline`` command, you can exit the session by pressing ``ctrl+b`` and then ``d``. This will detach the current tmux session from your terminal.

In order to reenter your tmux session, you can type in ``tmux ls`` in order to list all of your sessions and then ``tmux a -t <your-session-id>`` to re-enter that tmux session where your mitopipeline is running.