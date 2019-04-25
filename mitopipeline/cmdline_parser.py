import argparse, sys, os
from mitopipeline.pipeline_builder import PipelineBuilder
from mitopipeline.pipeline_runner import PipelineRunner
class CommandLineParser():

    def __init__(self, argv=sys.argv[1:]):
        self.__opts = self.parse_commands(argv)

    def parse_commands(self, argv=sys.argv[1:]):
        parser = self._build_parser()
        opts = parser.parse_args(argv)
        return opts

    def build_and_run(self, steps):
        pipeline_builder = PipelineBuilder()
        pipeline_builder.build_pipeline(slurm=self.__opts.slurm, tools=self.__opts.tools, directory=self.__opts.directory, steps=steps, output=self.__opts.output, refs=self.__opts.genomes)
        PipelineRunner.run(self.__opts)
        
    def _build_parser(self):
        parser = argparse.ArgumentParser()
        #required arguments
        required_args = parser.add_argument_group('required arguments')
        required_args.add_argument('-s', '--directory', help="Path to the directory of files to be run", type=str)
        required_args.add_argument('-g', '--genomes', help="Path of location of reference genomes", default=None)
        required_args.add_argument('-t', '--tools', help="Path to the directory that contains all of the 3rd party packages", default=None)
        #optional arguments
        parser.add_argument('-o', '--output', help="Path to where you want the output to be stored", default=None)
        parser.add_argument('-l', '--slurm', help="Use slurm jobs to run each step", default=False, action='store_true')
        parser.add_argument('-d', '--download', help="Specify softwares you want to download", default=False, action='store_true')
        parser.add_argument('-r', '--remove', nargs='+', help="Steps to not run in this pipeline", default=None)
        #parser.add_argument('-c', '--config', help="Use the config file to specify software options", default=None)
        parser.add_argument('-w', '--workers', type=int, help="Number of workers to use to run the pipeline", default=1)

        return parser
