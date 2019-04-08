import util
import luigi
import subprocess
import os
import shutil
STEPS = "./steps"
PIPELINE_START = "./"
PIPELINE_STORE = "./"

class PipelineRunner(luigi.WrapperTask):
    _file = luigi.Parameter(default="")
    def requires(self):
        if PIPELINE_START:
            for f in os.listdir(PIPELINE_START):
                if os.path.isfile(PIPELINE_START + f):
                    yield SNPEFF(id=util.parse_fid(f))
                    yield ANNOVAR(id=util.parse_fid(f))
                    yield HAPLOGREP(id=util.parse_fid(f))
                    
        else:
            yield SNPEFF(id=util.parse_fid(f))
            yield ANNOVAR(id=util.parse_fid(f))
            yield HAPLOGREP(id=util.parse_fid(f))
            
class SplitGap(luigi.Task):
    id=luigi.Parameter()
    def run(self):
        subprocess.call([STEPS + '/submit_splitgap_job.sh', self.id])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/SplitGap/' + '{}_split_gap'.format(self.id))

class RemoveNuMTs(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return SplitGap(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_removenumts_job.sh', self.id])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/RemoveNuMTs/' + '{}_remove_numts'.format(self.id))
class Downsample(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return RemoveNuMTs(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_downsample_job.sh', self.id])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/Downsample/' + '{}_downsample'.format(self.id))
class GATK(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return Downsample(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_gatk_job.sh', self.id])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/GATK/' + '{}_gatk'.format(self.id))
class SNPEFF(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return GATK(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_snpeff_job.sh', self.id])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/SNPEFF/' + '{}_snpeff'.format(self.id))
class ANNOVAR(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return GATK(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_annovar_job.sh', self.id])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/ANNOVAR/' + '{}_annovar'.format(self.id))
class HAPLOGREP(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return GATK(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_haplogrep_job.sh', self.id])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/HAPLOGREP/' + '{}_haplogrep'.format(self.id))
