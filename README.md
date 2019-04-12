### Tool Setup
This project can be downloaded through ``` pip install mitopipeline ``` followed by ``` pip install mitopipeline -r < requirements.txt ``` to install any dependencies. You can also directly download it through the source code at

Users will have to download 3rd party software packages that they wish to use on their own (such as GATK, ANNOVAR, etc.) before being able to use them in the pipeline.

Ensure that the directory where the tool is going to run on has the same naming convention. All files should be bam files and should end with a .bam extension. All steps within the pipeline will use the filename as the id, and add a _<stepname> to the filename id. For example, for filename ABCD.bam, the result of split_bam will be ABCD_splitbam.bam

### Documentation
Documentation can be viewed here: 

### Purpose
The purpose of this tool is to provide a useful pipeline for genetic researchers working with the mitochondrial genome. In genetics research, raw genome data requires many different softwares to be run first to extract useful information such as variant information before it can be analyzed. Thus, we want to provide a pipeline that allows users to streamline their processing automatically for them. 

### Users
This tool will be created for genetic researchers that have mitochondrial genome data in need of processing before analyzation through steps such as clipping, alignment into NuMT removals, running GATK, variant analysis, and other software packages such as SNPEFF, ANNOVAR, and Haplogrep.

### Goals
This tool should handle all dependency management in-between each step for the user from start to finish with the end result being processed data that the user can begin analyzing immediately. Specifically, this tool should:
- Handle all dep	endencies within the pipeline from start to finish
- Provide built in steps such as alignment, clipping, removal of NuMTs, and splitting the gap for Complete Genomics vcf files
- Use 3rd party software packages like GATK, ANNOVAR, SNPEFF, etc.
- Accurately log successes and failures that occurred while running the pipeline

### Acknowledgements
Many of the steps in the pipeline were adapted from scripts written by Sneha Grandhi. Check out her [Github](https://github.com/sneha-grandhi), [LinkedIn](https://www.linkedin.com/in/sneha-grandhi-phd-0165aa58/), or contact her through her [email](sneha_grandhi@hms.harvard.edu).
