Troubleshooting
***************

Outdated Modules on Case Western's HPC
--------------------------------------

As mentioned in the Running The Pipeline section, using mitopipeline on a separate server without admin privileges may require you to update various python modules on your own local ennvironment. For Case Western's HPC server, the details can be found `here https://sites.google.com/a/case.edu/hpcc/hpc-tutorials/installing-local-python-modules`_. In particular, python-dateutil will most likely need to be updated and installing mitopipeline itself will require you to use this method.

Reference Genomes
-----------------

Since the reference genomes are used by their exact file names, it's important to name all of your reference genomes the same way. Specifically, the reference genomes must be named hg38-nochr.fa, hg38.fa, and rCRS-MT.fa.

In addition, you must index the reference genomes through samtools/bwa. For each reference genome, the following complimentary file must exist (obtained through indexing):

<REF_name>.fa, <REF_name>.fa.amb, <REF_name>.fa.ann, <REF_name>.fa.bwt, <REF_name>.fa.fai, <REF_name>.fa.pac, <REF_name>.fa.sa, <REF_name>.dict

Lastly, the mitochondrial genome MUST be named "MT" (not chrM, chrM_rCRS, M) in all of the reference genomes. This is because in hg38, the mito genome is named MT, and matching to the reference contig is done through the literal string of the name.

SnpEff
------

SnpEff requires its own database and reference genome in order to be run. Thus, if it's your first time running the pipeline and you have multiple workers, it may be necessary to run 1 file first. Since SnpEff will be downloading its database the first time it runs, having multiple workers overlapping and trying to download the same database can cause some issues within the pipeline.

ANNOVAR
-------

Similar to SnpEff, ANNOVAR also requires its individual database and files in order to run. After registering and downloading, it's necessary to run table_anovar.pl before runnig mitopipeline. In particular, the `quick start http://annovar.openbioinformatics.org/en/latest/user-guide/startup/`_ page for ANNOVAR explains in detail how to do this. Note: the scripts in mitopipeline uses hg38, NOT hg19 as suggested on the website. Other than that, you can run each command as listed there.