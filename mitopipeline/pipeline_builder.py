import os, sys, ast, util, errors
from template_builder import task_template, task_with_req_template, import_template, paths_template, wrapper_task_template
SLURM = 1
NO_SLURM = 0
class PipelineBuilder():

    def __init__(self):
        self.softwares = ['gatk3', 'gatk4', 'snpeff', 'annovar', 'haplogrep']
        self.task_names = {'extract_mito': 'ExtractMito',
                            'split_gap': 'SplitGap',
                            'clipping': 'Clipping',
                            'remove_numts': 'RemoveNuMTs',
                            'downsample': 'Downsample',
                            'gatk': 'GATK',
                            'snpeff': 'SNPEFF',
                            'annovar': 'ANNOVAR',
                            'haplogrep': 'HAPLOGREP'}
        self.job_names = {'extract_mito': ['extractmito.sh', 'submit_extractmito_job.sh'],
                            'split_gap': ['splitgap.sh', 'submit_splitgap_job.sh'],
                            'clipping': ['clipping.sh', 'submit_clipping_job.sh'],
                            'remove_numts': ['remove_numts.sh', 'submit_removenumts_job.sh'],
                            'downsample': ['downsample.sh', 'submit_downsample_job.sh'],
                            'gatk3': ['gatk3.sh', 'submit_gatk3_job.sh'],
                            'gatk4': ['gatk4.sh', 'submit_gatk4_job.sh'],
                            'snpeff': ['snpeff.sh', 'submit_snpeff_job.sh'],
                            'annovar': ['annovar.sh', 'submit_annovar_job.sh'],
                            'haplogrep': ['haplogrep.sh', 'submit_haplogrep_job.sh']}

    def build_pipeline(self, tools=None, slurm=False, directory=None, steps=['extract_mito', 'split_gap', 'clipping', 'remove_numts', 'downsample', 'gatk', 'snpeff', 'annovar', 'haplogrep'], output=None):
        if not directory:
            raise ValueError('Building the pipeline requires a file/directory to run on')
        if not output:
            raise ValueError('Pipeline requires a directory to store output files')
        if not tools:
            raise ValueError('Pipeline builder requires a directory that contains all the software tools')
        else:
            errors.check_tools_exist(tools, steps)
            self.build_from_template(directory, output, steps, slurm)
    
    def build_from_template(self, directory, output, steps, slurm):
        errors.check_file_format(directory)
        errors.make_subdirectories(output, steps)
        with open('pipeline.py', 'w+') as pipeline:
            pipeline.write(import_template.render(imports=['util', 'luigi', 'subprocess', 'os', 'shutil']))
            pipeline.write(paths_template.render(directory=directory, output=output) + "\n\n")
            pipeline.write(wrapper_task_template.render(task_name="PipelineRunner", yields=list(self.task_names[step] for step in steps if step in self.softwares)) + "\n")
            prev_step = ""
            #if Complete Genomics data, i.e., did split gap then pipeline requires different scripts with shorter reads due to splitting into multiple reads at the gap
            if 'split_gap' in steps:
                self.job_names['remove_numts'] = ['remove_numts_no_split_gap.sh', 'submit_remove_numts_no_split_gap_job.sh']
                self.job_names['gatk3'] = ['gatk3_no_split_gap.sh', 'submit_gatk3_no_split_gap_job.sh']
                self.job_names['gatk4'] = ['gatk4_no_split_gap.sh', 'submit_gatk4_no_split_gap_job.sh']

            for step in steps:
                job_name = self.job_names[step][SLURM] if slurm else self.job_names[step][NO_SLURM]
                if prev_step == "":
                    pipeline.write(task_template.render(task_name=self.task_names[step], job_name=job_name, file_name=step) + "\n\n")
                    prev_step = self.task_names[step]
                else:
                    pipeline.write(task_with_req_template.render(task_name=self.task_names[step], req_name=prev_step, job_name=job_name, file_name=step) + "\n")
                    #softwares should only depend on gatk (the input file)
                    if "gatk" in step or step not in self.softwares:
                        prev_step = self.task_names[step]
                    
