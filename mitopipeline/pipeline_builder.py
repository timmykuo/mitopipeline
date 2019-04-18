import os, sys, pkg_resources
from mitopipeline.util import check_tools_exist, check_file_format, make_subdirectories, is_valid_directories, get_wrapper_tasks
from mitopipeline.templates import task_template, task_with_req_template, import_template, paths_template, wrapper_task_template, slurm_task_template, slurm_task_with_req_template

class PipelineBuilder():

    def __init__(self):
        self.FOLDER_NAME = 0
        self.EXTENSION = 1
        self.softwares = ['gatk', 'snpeff', 'annovar', 'haplogrep']
        self.task_template_info = {'extractmito': ['ExtractMito', 'extractmito.bam'],
                            'splitgap': ['SplitGap', 'splitgap.fastq'],
                            'clipping': ['Clipping', 'clipping.fastq'],
                            'removenumts': ['RemoveNuMTs', 'removenumts.bam'],
                            'downsample': ['Downsample', 'downsample.bam'],
                            'gatk': ['GATK', 'gatk.vcf'],
                            'snpeff': ['SNPEFF', 'snpeff.eff'],
                            'annovar': ['ANNOVAR', ''],
                            'haplogrep': ['HAPLOGREP', ''],
                            #for empty prevstep
                            '': ['', '']}
        self.dependencies = {'snpeff': ['snpeff'],
                        'annovar': ['annovar'],
                        'gatk': ['gatk', 'picard', 'rCRS'],
                        'removenumts': ['samtools', 'bwa', 'rCRS', 'hg38-nocrs'],
                        'clipping': ['samtools', 'bwa', 'bedtools', 'seqtk'],
                        'extractmito': ['samtools'],
                        'splitgap': ['samtools']}

    def build_pipeline(self, tools=None, slurm=False, directory=None, steps=['extractmito', 'splitgap', 'clipping', 'remove_numts', 'downsample', 'gatk', 'snpeff', 'annovar', 'haplogrep'], output=None, refs=None):
        is_valid_directories(directory, tools, refs, steps, self.softwares)
        check_tools_exist(tools, steps, self.dependencies)
        return self.build_from_template(directory, steps, slurm, output, tools, refs)

    def build_from_template(self, directory, steps, slurm, output, tools, refs):
        #user specified output or stored within mitopipeline directory
        output = output if output else "./pipeline_output"
        check_file_format(directory)
        make_subdirectories(output, self.task_template_info, steps, slurm)
        with open('./pipeline.py', 'w') as pipeline:
            pipeline.write(import_template.render(imports=['mitopipeline.util', 'luigi', 'subprocess', 'os', 'shutil', 'pkg_resources']))
            pipeline.write(paths_template.render(directory=directory, output=output, tools=tools, refs=refs) + "\n\n")
            pipeline.write(wrapper_task_template.render(task_name="PipelineRunner", yields=get_wrapper_tasks(self.task_template_info, steps, self.softwares)) + "\n")
            prev_step = ""
            #write in the steps requested into the pipeline
            for step in steps:
                #if Complete Genomics data, i.e., did split gap then pipeline requires different scripts with shorter reads due to splitting into multiple reads at the gap
                if 'splitgap' not in steps and step == 'removenumts':
                    job_name = 'remove_numts_no_split_gap.sh'
                #step only is name of the step, not the name of the script
                else:
                    job_name = step + ".sh"
                pipeline.write(self.get_template(slurm, prev_step, job_name, step))
                if "gatk" in step or step not in self.softwares:
                    prev_step = step
        return output
             
    def get_template(self, slurm, prev_step, job_name, step):
        task_name = self.task_template_info[step][self.FOLDER_NAME]
        file_name = self.task_template_info[step][self.EXTENSION]
        prev_file_name = self.task_template_info[prev_step][self.FOLDER_NAME]
        #if slurm job requested
        if slurm:
            #first step in the pipeline has no 'require' function
            if prev_step == "":
                return slurm_task_template.render(task_name=task_name, job_name=job_name, file_name=file_name) + "\n\n"
            #if not the first step in the pipeline, require the previous step
            else:
                return slurm_task_with_req_template.render(task_name=task_name, req_name=prev_file_name, job_name=job_name, file_name=file_name) + "\n"
        else:
            #first step in the pipeline has no 'require' function
            if prev_step == "":
                return task_template.render(task_name=task_name, job_name=job_name, file_name=file_name) + "\n\n"
            #if not the first step in the pipeline, require the previous step
            else:
                return task_with_req_template.render(task_name=task_name, req_name=prev_file_name, job_name=job_name, file_name=file_name) + "\n"
