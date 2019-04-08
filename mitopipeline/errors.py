import os, util
#checks that the file format follows our naming convenction
def check_file_format(directory):
    for f in os.listdir(directory):
        if not util.is_correct_format(f):
            raise ValueError(
                "All files saved in user-specified directory must follow the format 'NAME_ID_SAMPLETYPE'")

#check that all tools requested in steps are in the tools directory
def check_tools_exist(tools, steps, softwares):
    softwares = list(step for step in steps if step in softwares)
    for software in softwares:
        if not os.path.isdir(tools + "/" + software):
            raise ValueError(
                "User-specified 'tools' directory doesn't have a folder for " + software)

#creates subdirectories for all the requested steps within the specified output directory
def make_subdirectories(output, steps):
    for step in steps:
        if not os.path.isdir(output + "/" + step):
            os.makedirs(output + "/" + step)
