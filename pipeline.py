import mitopipeline.util
import luigi
import subprocess
import os
import shutil
import pkgutil
STEPS = os.path.join(os.path.dirname(__file__), 'steps')
PIPELINE_START = "../test/start"
PIPELINE_STORE = "../test/out"
SLURM_DIR = "../test/out/slurm"
TOOLS = "./tools"

class PipelineRunner(luigi.WrapperTask):
    def requires(self):
        for f in os.listdir(PIPELINE_START):
            if os.path.isfile(PIPELINE_START + "/" + f):
                yield RemoveNuMTs(id=mitopipeline.util.parse_fid(f))
                
class ExtractMito(luigi.Task):
    id=luigi.Parameter()
    def run(self):
        subprocess.call([STEPS + '/extract_mito.sh', self.id, PIPELINE_START, PIPELINE_STORE + '/ExtractMito', TOOLS])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/ExtractMito/' + '{}_extract_mito'.format(self.id))

class SplitGap(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return ExtractMito(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/split_gap.sh', self.id, PIPELINE_START + '/ExtractMito', PIPELINE_STORE + '/SplitGap', TOOLS])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/SplitGap/' + '{}_split_gap'.format(self.id))
class Clipping(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return SplitGap(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/clipping.sh', self.id, PIPELINE_START + '/SplitGap', PIPELINE_STORE + '/Clipping', TOOLS])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/Clipping/' + '{}_clipping'.format(self.id))
class RemoveNuMTs(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return Clipping(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/remove_numts_no_split_gap.sh', self.id, PIPELINE_START + '/Clipping', PIPELINE_STORE + '/RemoveNuMTs', TOOLS])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/RemoveNuMTs/' + '{}_remove_numts'.format(self.id))
