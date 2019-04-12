import os

def parse_fid(f):
    #if the file has been processed through at least one step in the pipeline
    if "_" in f:
        parsed = str(f).split("_")
    #filename is FILENAME.bam i.e., hasn't been processed yet
    else:
        parsed = str(f).split(".")
    return parsed[0]

def correct_format(f):
    return str(f).count('.') < 2

#checks that directory exists and mito directory contains the "steps" folder
def is_valid_directories(directory, tools, steps, softwares):
    if not directory:
        raise ValueError('Building the pipeline requires a file/directory to run on')
    check_tools_exist(tools, steps, softwares)
    #if not os.path.isdir(mito + "/steps"):
    #    raise ValueError('mito directory must map to the mitopipeline folder that contains /steps/')

#checks that the file format follows our naming convenction
def check_file_format(directory):
    for f in os.listdir(directory):
        if not correct_format(f):
            raise ValueError(
                "All files saved in user-specified directory must follow the format 'FILENAME.bam' with NO periods allowed in FILENAME")

#check that all tools requested in steps are in the tools directory
def check_tools_exist(tools, steps, softwares):
    softwares = list(step for step in steps if step in softwares)
    for software in softwares:
        if not tools:
            raise ValueError('Tools directory was not specified and is required if steps include a 3rd party software package')
        if not os.path.isdir(tools + "/" + software):
            raise ValueError(
                "User-specified 'tools' directory doesn't have a folder called " + software + " that contains the software")

#creates subdirectories for all the requested steps within the specified output directory
def make_subdirectories(output, steps, slurm):
    #create output folder that holds the mitopipeline output in the tool's directory
    if not os.path.isdir(output):
        os.makedirs(output)
    #TODO: fill in subdirectories for parts within each step
    subdirectories = {'remove_numts': ['fastqs', 'pileups', 'numt_removal_stor', 'counts'],
                        'split_gap': [],
                        'clipping': [],
                        'extract_mito': [],
                        'downsample': [],
                        'gatk': [],
                        'annovar': [],
                        'haplogrep': [],
                        'snpeff': [],
        }
    for step in steps:
        if not os.path.isdir(output + "/" + step):
            os.makedirs(output + "/" + step)
            for sub in subdirectories[step]:
                os.makedirs(output + "/" + sub)
    if slurm:
        os.makedirs(output + "/slurm")

def get_wrapper_tasks(task_names, steps, softwares):
    tasks = list(task_names[step] for step in steps if step in softwares)
    if not tasks:
        for task_name in reversed(list(task_names.keys())):
            #return the latest task that is not a software step
            if task_name not in softwares and task_name in steps:
                return [task_name]
    else:
        return tasks
