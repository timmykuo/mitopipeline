import os, sys
from mitopipeline.util import check_tools_exist, check_file_format, make_subdirectories, is_valid_directories
from mitopipeline.templates import task_template, task_with_req_template, import_template, paths_template, wrapper_task_template, slurm_task_template, slurm_task_with_req_template

class PipelineBuilder():

    def __init__(self):
        self.softwares = ['gatk', 'snpeff', 'annovar', 'haplogrep']
        self.task_names = {'extract_mito': 'ExtractMito',
                            'split_gap': 'SplitGap',
                            'clipping': 'Clipping',
                            'remove_numts': 'RemoveNuMTs',
                            'downsample': 'Downsample',
                            'gatk': 'GATK',
                            'snpeff': 'SNPEFF',
                            'annovar': 'ANNOVAR',
                            'haplogrep': 'HAPLOGREP'}

    def build_pipeline(self, tools=None, slurm=False, directory=None, steps=['extract_mito', 'split_gap', 'clipping', 'remove_numts', 'downsample', 'gatk', 'snpeff', 'annovar', 'haplogrep'], output=None):
        #TODO: placeholder for mito, need to figure out how to find the /steps folder when running from command line
        mito = "./"
        is_valid_directories(directory, tools, steps, self.softwares)
        self.build_from_template(directory, steps, slurm, mito, output, tools)

    def build_from_template(self, directory, steps, slurm, mito, output, tools):
        #user specified output or stored within mitopipeline directory
        output = output if output else mito + "/pipeline_output"
        check_file_format(directory)
        make_subdirectories(output, steps, slurm)
        with open('pipeline.py', 'w+') as pipeline:
            pipeline.write(import_template.render(imports=['mitopipeline.util', 'luigi', 'subprocess', 'os', 'shutil']))
            pipeline.write(paths_template.render(directory=directory, output=output, mito=mito, tools=tools) + "\n\n")
            pipeline.write(wrapper_task_template.render(task_name="PipelineRunner", yields=list(self.task_names[step] for step in steps if step in self.softwares)) + "\n")
            prev_step = ""
           
            #write in the steps requested into the pipeline
            for step in steps:
                #if Complete Genomics data, i.e., did split gap then pipeline requires different scripts with shorter reads due to splitting into multiple reads at the gap
                if 'split_gap' in steps and step == 'remove_numts':
                    job_name = 'remove_numts_no_split_gap.sh'
                #step only is name of the step, not the name of the script
                else:
                    job_name = step + ".sh"
                pipeline.write(self.get_template(slurm, prev_step, job_name, step))
                if "gatk" in step or step not in self.softwares:
                    prev_step = self.task_names[step]
             
    def get_template(self, slurm, prev_step, job_name, step):
        #if slurm job requested
        if slurm:
            #first step in the pipeline has no 'require' function
            if prev_step == "":
                return slurm_task_template.render(task_name=self.task_names[step], job_name=job_name, file_name=step) + "\n\n"
            #if not the first step in the pipeline, require the previous step
            else:
                return slurm_task_with_req_template.render(task_name = self.task_names[step], req_name = prev_step, job_name = job_name, file_name = step) + "\n"
        else:
            #first step in the pipeline has no 'require' function
            if prev_step == "":
                return task_template.render(task_name=self.task_names[step], job_name=job_name, file_name=step) + "\n\n"
            #if not the first step in the pipeline, require the previous step
            else:
                return task_with_req_template.render(task_name=self.task_names[step], req_name=prev_step, job_name=job_name, file_name=step) + "\n"
