import sys, os, subprocess

class PipelineRunner():

    @classmethod
    def run(cls, opts):
        if not os.path.isfile(opts.output + "/pipeline.py"):
            raise FileNotFoundError("Pipeline file was not found. Ensure that the directory and output paths are correct")
        else:
            command = "PYTHONPATH=\'" + str(opts.output) + "\' luigi --module pipeline PipelineRunner --workers " + str(opts.workers)+ " --local-scheduler"
            os.system(command)

    # @classmethod
    # def cleanup(cls, steps):