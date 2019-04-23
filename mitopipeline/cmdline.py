from mitopipeline.cmdline_parser import CommandLineParser
from mitopipeline.downloader import Downloader
#ALL_STEPS = ['extractmito', 'splitgap', 'clipping', 'removenumts', 'downsample',
#             'gatk', 'snpeff', 'annovar']
ALL_STEPS = ['extractmito', 'splitgap', 'clipping', 'removenumts', 'gatk', 'snpeff', 'annovar']

def run():
    cmdline_parser = CommandLineParser()
    opts = cmdline_parser.parse_commands()
    steps = remove_steps(opts.remove)
    if opts.download:
        downloader = Downloader()
        downloader.download_dependencies(steps)
    else:
        cmdline_parser.build_and_run(steps)

#Removes all steps listed in options
def remove_steps(steps):
    if not steps:
        return ALL_STEPS
    steps_to_use = ALL_STEPS
    for step in steps:
        if step in steps_to_use:
            steps_to_use.remove(step)
        else:
            raise ValueError(
                "The requested step to remove, " + step + ", doesn't exist in the pipeline")
    return steps_to_use
