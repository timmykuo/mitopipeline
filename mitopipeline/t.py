import luigi, subprocess
from os import path, listdir
from shutil import move
STEPS="./steps"
PIPELINE_START=""
PIPELINE_STORE=""

def parse_fid(f):
    parsed = str(f).split("_")
    cancer = parsed[0]
    patient_id = parsed[1] + "_" + parsed[2]
    return cancer + "_" + patient_id

def get_name(f_id):
    split = f_id.split("_")
    return split[0] + "_" + split[1]

# #run with options --workers to run luigi in parallel
# class RunDataPipeline(luigi.WrapperTask):
#     #parse_fid should return everything before .bam
#     #example filename: NBL_TARGET-30-PASMDM_01.bam
#     _file = luigi.Parameter(default="")

#     def requires(self):
#         if PIPELINE_START:
#             for f in listdir(PIPELINE_START):
#                 if path.isfile(PIPELINE_START + f):
#                     yield SNPEFF(id=parse_fid(f))
#                     yield ANNOVAR(id=parse_fid(f))
#                     yield HAPLOGREP(id=parse_fid(f))
#         else:
#             yield SNPEFF(id=parse_fid(f))
#             yield ANNOVAR(id=parse_fid(f))
#             yield HAPLOGREP(id=parse_fid(f))

# class ExtractMito(luigi.Task):
#     id = luigi.Parameter()

#     def run(self):
#         subprocess.call([STEPS + "submit_extractmito_job.sh", self.id])

#     def output(self):
#         return luigi.LocalTarget(PIPELINE_STORE + "/extracted_mito/" + "{}_mito.bam".format(self.id))
    
# class Clipping(luigi.Task):
#     id = luigi.Parameter()

#     def requires(self):
#         return ExtractMito(id=self.id)

#     def run(self):
#         subprocess.call([STEPS + "submit_clip_job.sh", self.id])

#     def output(self):
#         return luigi.LocalTarget(PIPELINE_STORE + "/clipped/" + "{}_clipped.bam".format(self.id))

# # This class splits all reads in a bam file wherever there is a gap based on the cigar string.
# #
# # Example: @ID TT with cigar string 1M2N1M would be split to the following reads:
# # @ID-A T 1M [other sam file columns]
# # @ID-B T 1M [other same file columns]
# #
# # Input: Bam file from /gdc-raw-mito-bams
# # Output: Fastq file that has been split
# class SplitGap(luigi.Task):
#     id = luigi.Parameter()

#     def requires(self):
#         return Clipping(id=self.id)

#     def run(self):
#         subprocess.call([STEPS + "submit_splitbam_job.sh", self.id])

#     def output(self):
#         return luigi.LocalTarget(PIPELINE_STORE + "/fastqs/" + "{}.fastq".format(self.id))

# class NuMTRemoval(luigi.Task):
#     id = luigi.Parameter()

#     def requires(self):
#         return SplitGap(id=self.id)

#     def run(self):
#         #currently only uses the updated fastq file
#         subprocess.call([STEPS + "submit_numtremoval_job.sh", self.id])

#     def output(self):
#         return luigi.LocalTarget(PIPELINE_STORE + "/numts/" + "{}.mito_hg38.sorted.mt.bam".format(self.id))


# class Downsample(luigi.Task):
#     id = luigi.Parameter()

#     def requires(self):
#         return NuMTRemoval(id=self.id)

#     def run(self):
#         subprocess.call([STEPS + "submit_downsample_job.sh", self.id])

#     def output(self):
#         return luigi.LocalTarget(PIPELINE_STORE + "/downsampled/" + "{}_downsampled.bam".format(self.id))

# class GATK(luigi.Task):
#     id = luigi.Parameter()

#     def requires(self):
#         return Downsample(id=self.id)

#     def run(self):
#         subprocess.call([STEPS + "submit_gatk_job.sh", self.id])

#     def output(self):
#         return luigi.LocalTarget(PIPELINE_STORE + "/gatk/" + "{}.gatk.vcf".format(self.id))


# class SNPEFF(luigi.Task):
#     id = luigi.Parameter()

#     def requires(self):
#         return GATK(id=self.id)

#     def run(self):
#         subprocess.call([STEPS + "submit_snpeff_job.sh", self.id])

#     def complete(self):
#         return luigi.LocalTarget(PIPELINE_STORE + "/snpeff/" + "{}_eff.vcf".format(self.id))

# class ANNOVAR(luigi.Task):
#     id = luigi.Parameter()

#     def requires(self):
#         return GATK(id=self.id)

#     def run(self):
#         subprocess.call([STEPS + "submit_annovar_job.sh", self.id])

#     def complete(self):
#         return luigi.LocalTarget(PIPELINE_STORE + "/annovar/" + "{}_annovar.vcf".format(self.id))


# class HAPLOGREP(luigi.Task):
#     id = luigi.Parameter()

#     def requires(self):
#         return GATK(id=self.id)

#     def run(self):
#         subprocess.call([STEPS + "submit_haplogrep_job.sh", self.id])

#     def complete(self):
#         return luigi.LocalTarget(PIPELINE_STORE + "/haplogrep/" + "{}_haplogrep.vcf".format(self.id))

