import argparse, sys, os
from pipeline_builder import PipelineBuilder
from pipeline_runner import PipelineRunner
class CommandLineParser():

    ALL_STEPS = ['extract_mito', 'split_gap', 'clipping', 'remove_numts', 'downsample',
                 'gatk', 'snpeff', 'annovar', 'haplogrep']
    DEPENDENCIES = {'snpeff': ['gatk'],
                'annovar': ['gatk'],
                'haplogrep': ['gatk'],
                'gatk': ['remove_numts'],
                'remove_numts': ['extract_mito', 'clip', 'split_gap']}

    def __init__(self, argv=sys.argv[1:]):
        self.__opts = self.parse_commands(argv)

    def parse_commands(self, argv=sys.argv[1:]):
        parser = self._build_parser()
        opts = parser.parse_args(argv)
        if not opts.directory:
            raise ValueError('Directory or file to be run on must be specified')
        return opts

    def build_and_run(self):
        pipeline_builder = PipelineBuilder()
        pipeline_builder.build_pipeline(slurm=self.__opts.slurm, tools=self.__opts.tools, directory=self.__opts.directory, steps=self.remove_steps(self.__opts.remove), output=self.__opts.output)
        #PipelineRunner.run('./pipeline.py', opts.workers)
        #if not opts.save:
        #    PipelineRunner.cleanup()
        
    def _build_parser(self):
        parser = argparse.ArgumentParser()
        #required arguments
        required_args = parser.add_argument_group('required arguments')
        required_args.add_argument('-d', '--directory', help="Path to the directory of files to be run", type=str)
        required_args.add_argument('-o', '--output', help="Output directory to save pipeline files")
        required_args.add_argument('-t', '--tools', help="Paths of locations of 3rd party packages")
        #optional arguments
        parser.add_argument('-s', '--save', help="Save files from the middle steps of pipeline instead of only the 3rd party software outputs", default=True, action='store_false')
        parser.add_argument('-l', '--slurm', help="Use slurm jobs to run each step", default=False, action='store_true')
        parser.add_argument('-r', '--remove', nargs='+', help="Steps to not run in this pipeline", default=None)
        #parser.add_argument('-c', '--config', help="Use the config file to specify software options", default=None)
        parser.add_argument('-w', '--workers', type=int, help="Number of workers to use to run the pipeline", default=1)

        return parser

    def remove_steps(self, steps):
        if not steps:
            return self.ALL_STEPS
        steps_to_use = self.ALL_STEPS
        for step in steps:
            if step in steps_to_use:
                steps_to_use.remove(step)
            else:
                raise ValueError("The requested step to remove doesn't exist in the pipeline")
        return steps_to_use
