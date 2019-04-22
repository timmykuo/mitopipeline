import os, shutil, pkg_resources
TOOLS = pkg_resources.resource_filename('mitopipeline', "tools")

def parse_fid(f):
    #filename is FILENAME.bam i.e.
    parsed = str(f).split(".")
    return parsed[0]

#ensure that all files in starting directory only has one period for the filename's extension, i.e. FILENAME.bam
def correct_format(f):
    return str(f).count('.') < 2

#checks that directory exists and mito directory contains the "steps" folder
def is_valid_directories(directory, tools, refs, steps, softwares):
    if not os.path.isdir(str(directory)):
        raise ValueError('Building the pipeline requires a file/directory to run on')
    if not os.path.isdir(str(refs)) and ("gatk" in steps or "removenumts" in steps):
        raise ValueError('GATK and RemoveNuMTs steps require a directory for the reference genomes')
    
#checks that the file format follows our naming convenction
def check_file_format(directory):
    for f in os.listdir(directory):
        #ignore hidden files
        if not f.startswith('.') and not correct_format(f):
            raise ValueError(
                "All files saved in user-specified directory must follow the format 'FILENAME.bam' with NO periods allowed in FILENAME")

#check that all tools required in steps are in the tools directory
def get_tools_loc(tools_dir, steps, dependencies):
    software_loc = {}
    softwares = set(dependencies[step] for step in steps)
    for software in softwares:
        loc = get_loc(software, tools_dir)
        software_loc[software] = loc
    return software_loc

def get_loc(software, tools_dir):
    #if available from command line
    if shutil.which(software):
        return 'command'
    #check if downloaded in mitopipeline's tools directory or subdirectory
    elif is_downloaded(software, TOOLS):
        return TOOLS
    elif is_downloaded(software, TOOLS + "/" + software):
        return TOOLS + "/" + software
    #check if downloaded in user-specified tools directory or subdirectory
    elif tools_dir:
        if is_downloaded(software, tools_dir):
            return tools_dir
        elif is_downloaded(software, tools_dir + "/" + software):
            return tools_dir + "/" + software
        else:
            raise ValueError("User-specified 'tools' directory doesn't have a folder called " + software + " that contains the software and that software is not available to run from the command line. Please install the required software through -d option or provide its' executable in " + tools_dir + "/" + software)
    else:
        raise ValueError("Software not available to run on command line. Please install the required software through -d option or provide a tools directory that contains its' executable in <path/to/tools_dir/" + software + "/>")


def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

#TODO: double check this function
def is_downloaded(program, tools):
    fpath, fname = os.path.split(tools)
    if fpath:
        if is_exe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return True

    return False

#creates subdirectories for all the requested steps within the specified output directory
def make_subdirectories(output, task_names, steps, slurm):
    #create output folder that holds the mitopipeline output in the tool's directory
    if not os.path.isdir(output):
        os.makedirs(output)
    #TODO: fill in subdirectories for parts within each step
    subdirectories = {'removenumts': ['fastqs', 'pileups', 'numt_removal_stor', 'counts'],
                        'splitgap': [],
                        'clipping': [],
                        'extractmito': [],
                        'downsample': [],
                        'gatk': ['gatk_stor'],
                        'annovar': [],
                        'haplogrep': [],
                        'snpeff': [],
        }
    for step in steps:
        folder_name = 0
        task_folder = output + "/" + task_names[step][folder_name]
        if not os.path.isdir(task_folder):
            os.makedirs(task_folder)
        for sub in subdirectories[step]:
            task_subfolder = task_folder + "/" + sub
            if not os.path.isdir(task_subfolder):
                os.makedirs(task_subfolder)
    if slurm:
        os.makedirs(output + "/slurm")

#returns either all of the softwares after gatk or the latest step before or
def get_wrapper_tasks(task_names, steps, softwares):
    folder_name = 0
    tasks = list(task_names[step][folder_name] for step in steps if step in softwares)
    if not tasks:
        for task_name in reversed(list(task_names.keys())):
            #return the latest task that is not a software step
            if task_name not in softwares and task_name in steps:
                #return the name of function in template instead of the step name
                return [task_names[task_name][folder_name]]
    else:
        return tasks
