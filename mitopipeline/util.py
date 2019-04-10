def parse_fid(f):
    #if the file has been processed through at least one step in the pipeline
    if "_" in f:
        parsed = str(f).split("_")
    #filename is FILENAME.bam i.e., hasn't been processed yet
    else:
        parsed = str(f).split(".")
    return parsed[0]
