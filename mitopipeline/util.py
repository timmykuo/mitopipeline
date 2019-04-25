import os, shutil, pkg_resources, subprocess, sys

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
def check_tools_exist(tools_dir, steps, dependencies):
    for step in steps:
        for dep in dependencies[step]:
            if not found_loc(dep, tools_dir):
                raise ValueError('Can\'t find ' + dep + ' in ' + tools_dir + ". Please download using -d option or make sure your tools directory has a folder called " + step)

def found_loc(software, tools_dir):
    #if available from command line
    if (software == 'samtools' or software == 'bwa'):
        if shutil.which(software):
            return True
        else:
            raise ValueError(software + " is not able to be run from the command line. Please refer to documentation on instructions for how to set up " + software + " or 'module load' it if your server uses Lmod")
    elif 'GenomeAnalysisTK' in software:
        return is_downloaded(software, tools_dir + "/gatk")
    elif 'snpEff' in software:
        return is_downloaded(software, tools_dir + "/snpEff")
    else:
        return is_downloaded(software, tools_dir + "/" + software)

def get_dir_name(software, dir):
    for name in os.listdir(dir):
        if os.path.isdir(dir + "/" + name) and software in name:
            return dir + "/" + name
    return None

def is_exe(fpath):
    return os.path.isfile(fpath)
    # and os.access(fpath, os.X_OK)

#checks if there is an executable called <program> on the path\
#software should be the software executable name
def is_downloaded(software, dir):
    print(dir + "/" + software)
    return is_exe(dir + "/" + software)

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
    #if snpeff and/or annovar are in tasks
    elif len(tasks) > 1 and "GATK" in tasks:
        tasks.remove("GATK")
    return tasks
        
def which(file):
    if shutil.which(file):
        return True
    for path in os.environ["PATH"].split(os.pathsep):
        if is_exe(os.path.join(path, file)):
            return True
    return False

def execute(cmd):
    popen = subprocess.Popen(cmd, stdout=subprocess.PIPE, universal_newlines=True)
    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line
    popen.stdout.close()
    return_code = popen.wait()
    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " [y/n] "
    elif default == "yes":
        prompt = " [Y/n] "
    elif default == "no":
        prompt = " [y/N] "
    else:
        raise ValueError("invalid default answer: '%s'" % default)

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

class cd:
    """Context manager for changing the current working directory"""

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)
