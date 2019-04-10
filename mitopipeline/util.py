import os

def parse_fid(f):
    #if the file has been processed through at least one step in the pipeline
    if "_" in f:
        parsed = str(f).split("_")
    #filename is FILENAME.bam i.e., hasn't been processed yet
    else:
        parsed = str(f).split(".")
    return parsed[0]

def is_correct_format(f):
    return str(f).count('.') >= 2

def is_valid_directories(directory, output, tools):
    if not directory:
        raise ValueError('Building the pipeline requires a file/directory to run on')
    if not output:
        raise ValueError('Pipeline requires a directory to store output files')
    if not tools:
        raise ValueError('Pipeline builder requires a directory that contains all the software tools')
    if not os.path.isdir("./steps"):
        raise ValueError('mitopipeline command must be run from the same folder that contains /steps/')

#checks that the file format follows our naming convenction
def check_file_format(directory):
    for f in os.listdir(directory):
        if not is_correct_format(f):
            raise ValueError(
                "All files saved in user-specified directory must follow the format 'FILENAME.bam'")

#check that all tools requested in steps are in the tools directory
def check_tools_exist(tools, steps, softwares):
    softwares = list(step for step in steps if step in softwares)
    for software in softwares:
        if not os.path.isdir(tools + "/" + software):
            raise ValueError(
                "User-specified 'tools' directory doesn't have a folder called " + software + " that contains the software")

#creates subdirectories for all the requested steps within the specified output directory
def make_subdirectories(output, steps, slurm):
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
