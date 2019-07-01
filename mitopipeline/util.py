import os, shutil, pkg_resources, subprocess, sys

#f is a file with an extension of .bam or .bai
def parse_fid(f):
    #remove extension
    parsed = str(f).split(".")
    fid = ""
    #ignore the extension
    for split in parsed[:-1]:
        fid += split
        fid += "."
    #remove the last period
    if fid != "":
        fid = fid[:-1]

    if f.endswith(("bam", "bai", "vcf")):
        return fid
    elif f.endswith("fastq"):
        if fid.endswith(("_1", "_2")):
            fid = fid[:-2]
        return fid
    else:
        raise ValueError("Parsing file must end with bam, bai, vcf, or fastq")

#ensure that all files in starting directory are same types of files
def correct_format(f, f_end):
    return f.endswith(f_end)

def get_f_end(f):
    if f.endswith("bam"):
        return "bam"
    elif f.endswith("bai"):
        return "bai"
    elif f.endswith("fastq"):
        return "fastq"
    else:
        return ValueError('Files must be either bam or fastq files')

#checks that directory exists and mito directory contains the "steps" folder
def is_valid_directories(directory, refs, steps):
    if not os.path.isdir(str(directory)):
        raise ValueError('Building the pipeline requires a file/directory to run on')
    if not os.path.isdir(str(refs)) and ("gatk" in steps or "removenumts" in steps):
        raise ValueError('GATK and RemoveNuMTs steps require a directory for the reference genomes')
    else:
        return True
    
def correct_start_files(step, directory):
    if (step == "extractmito" or step == "clipping" or step == "splitgap" or step == "gatk") and not all_same_end(directory, ("bam", "bai")):
        raise ValueError("First step in your pipeline is " + step + ". All files in your start directory must end with bam/bai")
    elif (step == "removenumts") and not all_same_end(directory, ("fastq")):
        raise ValueError("First step in your pipeline is " + step + ". All files in your start directory must end with fastq")
    elif (step == "snpeff" or step == "annovar") and not all_same_end(directory, ("vcf")):
        raise ValueError("First step in your pipeline is " + step + ". All files in your start directory must end with vcf")
    return True

def all_same_end(directory, ends):
    for f in os.listdir(directory):
        if not f.endswith(ends):
            return False
    return True

#checks that the file format follows our naming convenction
def files_correct_format(directory):
    for f in os.listdir(directory):
        f_end = get_f_end(f)
        if not f.startswith('.') and not correct_format(f, f_end):
            raise ValueError(
                "All files saved in user-specified directory must follow the same format of 'FILENAME.bam' or 'FILENAME_1.fastq/FILENAME_2.fastq/FILENAME.fastq' (depending on if it's paired end).")
    return True

#check that all tools required in steps are in the tools directory
def tools_exist(tools_dir, steps, dependencies):
    for step in steps:
        for dep in dependencies[step]:
            if not found_loc(dep, tools_dir):
                raise ValueError('Can\'t find ' + dep + ' in ' + tools_dir + ". Please download using -d option or make sure your tools directory has a folder called " + dep)
    return True

#function to check annovar dependencies, can't make it general since the files are specific
def is_annovar_downloaded(software, tools_dir):
    return os.path.isfile(tools_dir + "/" + software)
    
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
    elif 'annovar' in software:
        return is_annovar_downloaded(software, tools_dir + "/annovar")
    else:
        return is_downloaded(software, tools_dir + "/" + software)

def get_dir_name(software, dir):
    for name in os.listdir(dir):
        if os.path.isdir(dir + "/" + name) and software in name:
            return dir + "/" + name
    return None

def is_exe(fpath):
    return os.path.isfile(fpath) and os.access(fpath, os.R_OK)

#checks if there is an executable called <program> on the path\
#software should be the software executable name
def is_downloaded(software, dir):
    return is_exe(dir + "/" + software)

#creates subdirectories for all the requested steps within the specified output directory
def make_subdirectories(output, task_names, steps, email):
    #create output folder that holds the mitopipeline output in the tool's directory
    if not os.path.isdir(output):
        os.makedirs(output)
    subdirectories = {'removenumts': ['fastqs', 'pileups', 'numt_removal_stor', 'counts'],
                        'gatk': ['gatk_stor']}
    for step in steps:
        folder_name = 0
        task_folder = output + "/" + task_names[step][folder_name]
        if not os.path.isdir(task_folder):
            os.makedirs(task_folder)
        if step in subdirectories:
            for sub in subdirectories[step]:
                task_subfolder = task_folder + "/" + sub
                if not os.path.isdir(task_subfolder):
                    os.makedirs(task_subfolder)
    if email and not os.path.isdir(output + "/slurm"):
        os.makedirs(output + "/slurm")
        os.makedirs(output + "/slurm/STDOUT")

#returns either all of the softwares after gatk or the latest step before or
def get_wrapper_tasks(task_names, steps, softwares):
    folder_name = 0
    #tasks is a list of the steps' foldernames included in softwares
    tasks = list(task_names[step][folder_name] for step in steps if step in softwares)
    if not tasks:
        #find the latest task that is not a software step
        for task_name in reversed(list(task_names.keys())):
            if task_name not in softwares and task_name in steps:
                #return foldername of the task
                return [task_names[task_name][folder_name]]
    if "GATK" not in tasks:
        raise ValueError("GATK must be in the list of steps if snpeff and/or annovar are included")
    #snpeff and/or annovar are in tasks, yield those instead of gatk
    if len(tasks) > 1:
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
