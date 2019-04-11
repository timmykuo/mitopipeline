import mitopipeline.util
import luigi
import subprocess
import os
import shutil
STEPS = ".//steps"
PIPELINE_START = "../../test/start"
PIPELINE_STORE = ".//pipeline_output"
SLURM_DIR = ".//pipeline_output/slurm"
TOOLS = "./tools"

class PipelineRunner(luigi.WrapperTask):
    _file = luigi.Parameter(default="")
    def requires(self):
        if PIPELINE_START:
            for f in os.listdir(PIPELINE_START):
                if os.path.isfile(PIPELINE_START + f):
                    yield GATK(id=mitopipeline.util.parse_fid(f))
                    
        else:
            yield GATK(id=mitopipeline.util.parse_fid(f))
            
class SplitGap(luigi.Task):
    id=luigi.Parameter()
    def run(self):
        subprocess.call([STEPS + '/split_gap.sh', self.id, PIPELINE_START, PIPELINE_STORE + '/SplitGap', TOOLS])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/SplitGap/' + '{}_split_gap'.format(self.id))

class RemoveNuMTs(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return SplitGap(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/remove_numts_no_split_gap.sh', self.id, PIPELINE_START + '/SplitGap', PIPELINE_STORE + '/RemoveNuMTs', TOOLS])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/RemoveNuMTs/' + '{}_remove_numts'.format(self.id))
class GATK(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return RemoveNuMTs(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/gatk.sh', self.id, PIPELINE_START + '/RemoveNuMTs', PIPELINE_STORE + '/GATK', TOOLS])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/GATK/' + '{}_gatk'.format(self.id))
