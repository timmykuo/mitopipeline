import util
import luigi
import subprocess
import os
import shutil
STEPS = "./steps"
PIPELINE_START = "./"
PIPELINE_STORE = "./"
SLURM_DIR = "./slurm"

class PipelineRunner(luigi.WrapperTask):
    _file = luigi.Parameter(default="")
    def requires(self):
        if PIPELINE_START:
            for f in os.listdir(PIPELINE_START):
                if os.path.isfile(PIPELINE_START + f):
                    yield GATK(id=util.parse_fid(f))
                    yield SNPEFF(id=util.parse_fid(f))
                    yield ANNOVAR(id=util.parse_fid(f))
                    yield HAPLOGREP(id=util.parse_fid(f))
                    
        else:
            yield GATK(id=util.parse_fid(f))
            yield SNPEFF(id=util.parse_fid(f))
            yield ANNOVAR(id=util.parse_fid(f))
            yield HAPLOGREP(id=util.parse_fid(f))
            
class ExtractMito(luigi.Task):
    id=luigi.Parameter()
    def run(self):
        subprocess.call([STEPS + '/submit_job.sh', self.id, extract_mito.sh, SLURM_DIR, PIPELINE_START, PIPELINE_STORE + '/ExtractMito'])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/ExtractMito/' + '{}_extract_mito'.format(self.id))

class SplitGap(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return ExtractMito(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_job.sh', self.id, split_gap.sh, SLURM_DIR, PIPELINE_START + '/ExtractMito', PIPELINE_STORE + '/SplitGap'])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/SplitGap/' + '{}_split_gap'.format(self.id))
class Clipping(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return SplitGap(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_job.sh', self.id, clipping.sh, SLURM_DIR, PIPELINE_START + '/SplitGap', PIPELINE_STORE + '/Clipping'])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/Clipping/' + '{}_clipping'.format(self.id))
class RemoveNuMTs(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return Clipping(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_job.sh', self.id, remove_numts.sh, SLURM_DIR, PIPELINE_START + '/Clipping', PIPELINE_STORE + '/RemoveNuMTs'])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/RemoveNuMTs/' + '{}_remove_numts'.format(self.id))
class Downsample(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return RemoveNuMTs(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_job.sh', self.id, downsample.sh, SLURM_DIR, PIPELINE_START + '/RemoveNuMTs', PIPELINE_STORE + '/Downsample'])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/Downsample/' + '{}_downsample'.format(self.id))
class GATK(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return Downsample(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_job.sh', self.id, gatk.sh, SLURM_DIR, PIPELINE_START + '/Downsample', PIPELINE_STORE + '/GATK'])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/GATK/' + '{}_gatk'.format(self.id))
class SNPEFF(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return GATK(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_job.sh', self.id, snpeff.sh, SLURM_DIR, PIPELINE_START + '/GATK', PIPELINE_STORE + '/SNPEFF'])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/SNPEFF/' + '{}_snpeff'.format(self.id))
class ANNOVAR(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return GATK(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_job.sh', self.id, annovar.sh, SLURM_DIR, PIPELINE_START + '/GATK', PIPELINE_STORE + '/ANNOVAR'])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/ANNOVAR/' + '{}_annovar'.format(self.id))
class HAPLOGREP(luigi.Task):
    id=luigi.Parameter()
    def require(self):
        return GATK(id=self.id)
    def run(self):
        subprocess.call([STEPS + '/submit_job.sh', self.id, haplogrep.sh, SLURM_DIR, PIPELINE_START + '/GATK', PIPELINE_STORE + '/HAPLOGREP'])
    def output(self):
        return luigi.LocalTarget(PIPELINE_STORE + '/HAPLOGREP/' + '{}_haplogrep'.format(self.id))
