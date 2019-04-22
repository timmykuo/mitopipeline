Example - One TARGET Genome File
********************************

This section will go through an example of how to run the pipeline step by step assuming that you've gone ahead and installed it through pip.

Step1 - Download steps if necessary
-----------------------------------

Step2 - Find subdirectories such as starting directory, out directory, REFs and Tools
-------------------------------------------------------------------------------------

In the example case, it's located within /mitopipeline/example, so depending on where your mitopipeline package is installed, you will have to specify that directory. 

In addition, the REFs directory will need to be referenced based on where you have the human reference genomes. The requirements for when the reference genomes are needed can be found in the pipeline steps page. In this case, the tools directory specified will just be in the same folder as the start and out folders. Thus, the command I will be running is:

Thus, our following directories to be specified are:

- start directory == './example'
- out directory == '../test/out'
- tools directory == './tools'
- refs directory == './REFs'

Step3 - Run pipeline
--------------------

Now that we have our directories, we can run the pipeline. 

.. code:: bash

    $ mitopipeline -s ./example/ -o ./test/out -g ./REFs/ -t ./tools -r clipping extractmito

One thing to note is that since this example file is from the TARGET database and was sequenced by Complete Genomics, we will be using the split gap step but without clipping since our reads are relatively short already. Also, since our starting file is alreay a mitochondrial genome, we can also skip the extract mito step.

It's also worth noting that we do not specify the number of workers. Since we only have one file in our example folder, it's not necessary to have multiple workers and so we can use the default of 1.

Step4 - View Output
-------------------

Once the pipeline has finished running, you will see either a successful or error. You can read the log output to see what parts may have failed and which have succeeded. All output files can be viewed within the specified out directory with subdirectories created for each step.