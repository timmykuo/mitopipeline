import sys, os

class PipelineRunner():

    @classmethod
    def run(cls, pipeline_loc, opts):
        if not os.path.isfile(pipeline_loc + "/pipeline.py"):
            raise FileNotFoundError("Pipeline file was not found. Ensure that the directory and output paths are correct")
        else:
            command = "PYTHONPATH=\'"+pipeline_loc+"\' luigi --module pipeline PipelineRunner --workers " + str(opts.workers)+ " --local-scheduler"
            os.system(command)

    # @classmethod
    # def cleanup(cls, steps):



