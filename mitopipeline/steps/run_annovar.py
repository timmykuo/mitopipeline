import os, subprocess

def parse_fid(f):
    split = f.split(".")
    return split[0]

def is_empty(f):
    with open("/mnt/rds/txl80/LaframboiseLab/tyk3/scripts/gatk3_vcfs/" + f) as vcf:
        size = os.path.getsize("/mnt/rds/txl80/LaframboiseLab/tyk3/scripts/gatk3_vcfs/" + f)
        for line in vcf:
            size -= len(line)
            if not size and "CHROM" in line and "POS" in line and "ID" in line and "REF" in line and "ALT" in line and "QUAL" in line and "FILTER" in line and "FORMAT" in line and "foo" in line:
                return True
    return False

def run_annovar():
    for f in os.listdir("/mnt/rds/txl80/LaframboiseLab/tyk3/scripts/gatk3_vcfs"):
        fid = parse_fid(f)
        if not is_empty(f):
            subprocess.call(["/mnt/rds/txl80/LaframboiseLab/tyk3/scripts/run_annovar.sh", fid])

if __name__ == "__main__":
    run_annovar()
