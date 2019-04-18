Pipeline Steps
****************

A number of steps are able to be included within the pipeline that can be controlled through options in the command line. As mentioned on the home page, many of these files were adapted from bash scripts written by Sneha Grandhi. 

All of them utilize samtools, bwa or both and thus are necessary to download. This section will explain each step in further detail so you can determine if it's necessary to include within your pipeline.

Submit Job
----------

Extract Mito
------------

Requirements: samtools

.. code:: bash

    samtools index $START/$FILENAME.bam
    samtools view -H $START/$FILENAME.bam > $START/$FILENAME.header.bam
    grep -i SN:MT $START/$FILENAME.header.bam > $START/$FILENAME.MT.bam
    grep -i SN:chrM_rCRS $START/$FILENAME.header.bam > $START/$FILENAME.chrM_rCRS.bam
    grep -i SN:chrM $START/$FILENAME.header.bam > $START/$FILENAME.chrM.bam
    grep -i SN:M $START/$FILENAME.header.bam > $START/$FILENAME.M.bam

This snippet extracts the mitochondrial genome from the bam file that was passed in where ``$START`` is the directory that contains the bam files and ``$FILENAME`` is the filename. Since the mitochondrial genome can have different notations based on which reference genome it was aligned to, the script takes into account all the different mitochondrial genome tags.

Clipping
--------

Requirements: samtools, bwa, bam2fastq, seqtk-master



SplitGap
--------

Requirements: samtools 

Remove NuMTs
------------

Requirements: samtools, bwa, bam2fastq, hg38 mitochondrial reference genome (rCRS-MT.fa), hg38 human genome without mitochondrial genome (hg38-norcrs.fa)

GATK
----

Requirements: GATK-3.1, picard-1.93, hg38 mitochondrial reference genome (rCRS-MT.fa)

SNPEFF
------

Requirements: snpEff

ANNOVAR
-------

Requirements: Annovar



