import configparser

args = ['directory', 'file', 'output', 'steps', 'workers']

#required arguments
files = {'directory': None, 'file': None}
output = None

#optional arguments
steps = ['extract_mito', 'split_gap', 'clip', 'remove_numts', 'downsample',
         'gatk', 'heteroplasmy_analysis', 'gatk', 'snpeff', 'annovar', 'haplogrep']
workers = 1

def set_arg(arg, val):
    if arg not in args:
        raise ValueError('Input argument is not an available option')
    elif arg == "directory" or arg == "file":
        files.arg = val
    else:
        files.arg = val

