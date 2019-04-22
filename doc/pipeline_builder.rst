Pipeline Steps
****************

A number of steps are able to be included within the pipeline that can be controlled through options in the command line. As mentioned on the home page, many of these files were adapted from bash scripts written by Sneha Grandhi. 

All of them utilize samtools, bwa or both and thus are necessary to download. This section will explain each step in further detail so you can determine if it's necessary to include within your pipeline. 

If you decide to use the `-l` (slurm) option, then each file that is run on each step of the pipeline will be submitted as a slurm job. Options that you will need to specify are described in the command line options section.

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

Requirements: samtools, bwa, bam2fastq, seqtk

The clipping step takes in as input a bam file, either from the mitochondrial bam file extracted from the previous step, or as a starting point. 

First, it changes the bam file into a fastq format:

.. code:: bash

    $TOOLS/bam2fastq -f -o $OUT/$1_bam2fastq_#.fastq $START/$1_$filetype.bam

Then, depending on if it's a paired end or single end, will trim the first and last two pairs from every read.

.. code:: bash

    if [ -e "$OUT/$1_bam2fastq_1.fastq" ]
    then
    echo "PAIRED-END"
            echo "--CLIPPED: Removing first and last 2 base pairs from every read"
            $TOOLS/seqtk-master/seqtk trimfq -b 2 -e 2 $OUT/$1_bam2fastq_1.fastq > $OUT/$1_1.fastq
            $TOOLS/seqtk-master/seqtk trimfq -b 2 -e 2 $OUT/$1_bam2fastq_2.fastq > $OUT/$1_2.fastq
    else
    echo "SINGLE-END"
            echo "--CLIPPED: Removing first and last 2 base pairs from every read"
            $TOOLS/seqtk-master/seqtk trimfq -b 2 -e 2 $OUT/$1_bam2fastq.fastq > $OUT/$1.fastq
    fi

The output of this step will be $filename.fastq or $filename_1.fastq and $filename_2.fastq, depending on if it's a paired end read or not.

SplitGap
--------

Requirements: samtools 

This step is particularly for Complete Genomics data. Complete Genomics uses a technique called `DNA Nanoball Sequencing <https://en.wikipedia.org/wiki/DNA_nanoball_sequencing>`_. 

The downside to this technique is that each read ends up having a 'gap' in the middle of the read that doesn't match to the reference genome. For example, a reference genome read of 'ATCGA' on a DNA Nanoball Sequencing read could be 'ATA' where the cigar string would be 2M2N1M. To fix this, we use the sam file's format `cigar string <https://www.drive5.com/usearch/manual/cigar.html/>`_ to split the read where the matches are into two reads. So the previous examples reads would be split into 'AT' and 'A' with 2M and 1M being their respective cigar strings.

.. code:: python

    #regex matches any amount of numbers followed by one capital letter
    matches = re.findall(r'[0-9]+[A-Z]{1}', cigar)
    curr_num = 0
    temp_idx = 0
    let = ""
    for match in matches:
        num = int(match[:-1])
        let = match[-1]
        #append next reads and quals if not an M or an I
        if let != 'M' or let != 'I' and curr_num != 0:
            add_next_seq(temp_idx, curr_num, new_reads, reads, new_quals, quals)
            temp_idx += curr_num
            curr_num = 0
        else:
            curr_num += num

The above code block finds all cigar blocks, i.e. any number followed by a capital letter like 1M, 5N, 40I, etc.  It then appends on the next reads only if the letter is not an M or an I, which are the reads that match up to the reference genome.

Although this does decrease the read size since they are being split, the quality of the reads drastically improve.

Remove NuMTs
------------

Requirements: samtools, bwa, bam2fastq, hg38 mitochondrial reference genome (rCRS-MT.fa), hg38 human genome without mitochondrial genome (hg38-norcrs.fa), and hg38 human genome (hg38.fa)

`NuMTs <https://en.wikipedia.org/wiki/NUMT>`_ are DNA sequences harbored in the nuclear genome, but closely resemble sequences in the mitochondrial genome. We remove these as quality control and to reduce noise in the following steps. The output of this step is a bam file with NuMTs removed

To do this, we first align our input fastq files to both the mitochondrial genome and hg38 without the mitochondrial genome to find any close matches. Then, we extract the perfect matches to the nuclear genome, realign the resulting fastq file back to hg38 reference genome, and extract the mitochondrial genome. 

GATK
----

Requirements: gatk-3.1, picard, hg38 mitochondrial reference genome (rCRS-MT.fa)

The gatk script were adapted from the suggested pipeline by GATK. In particular, the following steps are run in order:

Picard's AddOrReplaceReadGroups, Picard's MarkDuplicates, GATK's RealignerTargetCreator, GATK's IndelRealigner, GATK's FixMateInformation, GATK's BaseRecalibrator, GATK's PrintReads, GATK's HaplotypeCaller, GATK's VariantFiltration.

An example of how gatk is called:

.. code:: bash

    java -Xmx10g -jar $TOOLS/gatk/gatk.jar \
    -T HaplotypeCaller \
    -R $REFS/rCRS-MT.fa \
    -I $TMPDIR/$1.tcga.marked.realigned.fixed.read.bam \
    --maxReadsInRegionPerSample 200 \
    --sample_ploidy 100 \
    -stand_call_conf 50 \
    -stand_emit_conf 10 \
    --pcr_indel_model HOSTILE \
    -minPruning 10 \
    -A StrandAlleleCountsBySample \
    --dbsnp $5/dbsnp/mtONLY.vcf \
    -o $TMPDIR/$1.tcga.snps.vcf

Something important to note is that the gatk.jar executable must be placed within a folder called gatk within the tool's directory.

SNPEFF
------

Requirements: snpEff

This use's snpeff's most basic command and using the most recent mitochondrial reference genome GRCh38.86

.. code:: bash

    java -Xmx4g -jar $TOOLS/snpEff/snpEff.jar GRCh38.86 $VCFS/$1_$filetype.vcf > $SNPEFF/$1_snpEff.vcf

This is the standard usage of snpEff. You can read more about it on their website. Also note that the snpEff executable must be placed within a snpEff folder within the tool's directory just like gatk.

ANNOVAR
-------

Requirements: Annovar

Annovar can only be downloaded after registering on their `website <http://www.openbioinformatics.org/annovar/annovar_download_form.php>`_.

.. code:: bash

    #convert vcf file to avinput file
    perl $TOOLS/convert2annovar.pl -format vcf4 $VCFS/$1_$filetype.vcf  > $ANNOVAR/$1.avinput

    perl $TOOLS/table_annovar.pl $ANNOVAR/$1.avinput $TOOLS/humandb/ -remove -protocol dbnsfp33a -operation f -build hg38 -nastring . > $3/$1.avoutput

This is the suggested usage from annovar. You can read more about these files on their website. It also uses the hg38 version of the human reference genome.




