import os, sys, pkg_resources
from mitopipeline.util import tools_exist, files_correct_format, make_subdirectories, is_valid_directories, get_wrapper_tasks, correct_start_files, all_same_end
from mitopipeline.templates import task_template, task_with_req_template, import_template, paths_template, wrapper_task_template, slurm_task_template, slurm_task_with_req_template

class PipelineBuilder():

    def __init__(self):
        #hack to find tools folder in cache
        try:
            self.TOOLS = pkg_resources.resource_filename('mitopipeline', "tools")
        except KeyError:
            self.TOOLS = pkg_resources.resource_filename('mitopipeline', "/") + "/tools"
        if not os.path.isdir(self.TOOLS):
            os.makedirs(self.TOOLS)
        self.FOLDER_NAME = 0
        self.EXTENSION = 1
        self.softwares = ['gatk', 'snpeff', 'annovar']
        #folder output of the step, and output extension added onto the fild id of each step, ex. extractmito step --> filename_extractmito.bam
        self.task_template_info = {'extractmito': ['ExtractMito', 'extractmito.bam'],
                            'splitgap': ['SplitGap', '1_splitgap.fastq'],
                            'clipping': ['Clipping', '1_clipping.fastq'],
                            'removenumts': ['RemoveNuMTs', 'removenumts.bam'],
                            'downsample': ['Downsample', 'downsample.bam'],
                            'gatk': ['GATK', 'gatk.vcf'],
                            'snpeff': ['SNPEFF', 'snpeff.vcf'],
                            'annovar': ['ANNOVAR', 'avinput.hg38_multianno.txt'],
                            #for empty prevstep
                            '': ['', '']}
        #software dependencis for each step of the pipeline
        self.dependencies = {'snpeff': ['snpEff.jar'],
                             'gatk': ['GenomeAnalysisTK.jar'],
                             'annovar': ['table_annovar.pl', 'convert2annovar.pl'],
                             'removenumts': ['samtools', 'bwa'],
                             'clipping': ['samtools', 'bwa', 'seqtk'],
                             'extractmito': ['samtools'],
                             'splitgap': ['samtools']}

    def build_pipeline(self, tools=None, directory=None, steps=['extractmito', 'splitgap', 'clipping', 'remove_numts', 'gatk', 'snpeff', 'annovar'], output=None, refs=None, email=None):
        tools = tools if tools else self.TOOLS
        if is_valid_directories(directory, refs, steps) and tools_exist(tools, steps, self.dependencies) and correct_start_files(steps[0], directory):
            return self.build_from_template(directory, steps, output, tools, refs, email)
        else:
            raise ValueError("Errors in directories and/or software requirements.")

    def build_from_template(self, directory, steps, output, tools, refs, email):
        #user specified output or stored within mitopipeline directory
        output = output if output else "./pipeline_output"
        if files_correct_format(directory):
            make_subdirectories(output, self.task_template_info, steps, email)
            with open(output + '/pipeline.py', 'w') as pipeline:
                pipeline.write(import_template.render(imports=['mitopipeline.util', 'luigi', 'subprocess', 'os', 'shutil', 'pkg_resources']))
                pipeline.write(paths_template.render(directory=directory, output=output, tools=tools, refs=refs) + "\n\n")
                pipeline.write(wrapper_task_template.render(task_name="PipelineRunner", yields=get_wrapper_tasks(self.task_template_info, steps, self.softwares)) + "\n")
                prev_step = ""
                #write in the steps requested into the pipeline
                for step in steps:
                    #step only is name of the step, not the name of the script
                    job_name = step + ".sh"
                    #if Complete Genomics data, i.e., did split gap then pipeline requires different scripts with shorter reads due to splitting into multiple reads at the gap
                    if 'splitgap' not in steps:
                        if step == 'removenumts':
                            job_name = 'remove_numts_no_split_gap.sh'
                    pipeline.write(self.get_template(prev_step, job_name, step, email))
                    if "gatk" in step or step not in self.softwares:
                        prev_step = step
        return os.path.isfile(output + "/pipeline.py")
             
    def get_template(self, prev_step, job_name, step, email):
        task_name = self.task_template_info[step][self.FOLDER_NAME]
        file_name = self.task_template_info[step][self.EXTENSION]
        prev_task_name = self.task_template_info[prev_step][self.FOLDER_NAME]
        #if slurm job requested
        if email:
            #first step in the pipeline has no 'require' function
            if prev_step == "":
                return slurm_task_template.render(task_name=task_name, job_name=job_name, file_name=file_name, email=email) + "\n\n"
            #if not the first step in the pipeline, require the previous step
            else:
                return slurm_task_with_req_template.render(task_name=task_name, req_name=prev_task_name, job_name=job_name, file_name=file_name, email=email) + "\n"
        else:
            #first step in the pipeline has no 'require' function
            if prev_step == "":
                return task_template.render(task_name=task_name, job_name=job_name, file_name=file_name) + "\n\n"
            #if not the first step in the pipeline, require the previous step
            else:
                return task_with_req_template.render(task_name=task_name, req_name=prev_task_name, job_name=job_name, file_name=file_name) + "\n"
