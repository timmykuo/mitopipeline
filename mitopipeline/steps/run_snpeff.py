import os, subprocess, time
SNPEFFS = "/mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/snpeffs"
#VCFS = "/mnt/rds/txl80/LaframboiseLab/tyk3/TARGET_PIPELINE/gatk_vcfs"
VCFS = "/mnt/rds/txl80/LaframboiseLab/tyk3/scripts/gatk3_vcfs"
def parse_fid(f):
    split = f.split(".")
    return split[0]

def run_snpeff():
    for f in os.listdir(VCFS):
        f_id = parse_fid(f)
        f_name = f_id + "_snpEff.vcf"
        if not f_name in os.listdir(SNPEFFS):
            subprocess.call(["./submit_snpeff_job.sh", f_id])

if __name__ == "__main__":
    while True:
        run_snpeff()
