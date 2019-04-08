def parse_fid(f):
    parsed = str(f).split("_")
    cancer = parsed[0]
    patient_id = parsed[1] + "_" + parsed[2]
    return cancer + "_" + patient_id


def get_name(f_id):
    split = f_id.split("_")
    return split[0] + "_" + split[1]
