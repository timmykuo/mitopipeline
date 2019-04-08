import sys, os

class PipelineRunner():

    @classmethod
    def run(cls, filename, opts):
        if not os.path.isfile(filename):
            raise FileNotFoundError("Pipeline file was not found. Ensure that the directory and output paths are correct")
        else:
            command = "PYTHONPATH=\'.\' luigi --module pipeline RunDataPipeline --workers " + opts.workers + " --local-scheduler"
            os.system(command)



