import sys
def get_read_depths(fastq_file):
    with open(fastq_file) as fastq:
        depth = 0
        num_lines = 0
        line_num = 1
        for line in fastq:
            #2nd line, then get every 4th line
            if line_num == 2 or (line_num - 2) % 4 == 0:
                depth += len(line)
                num_lines += 1
            line_num += 1
        if num_lines == 0:
            avg_depth = 0
        else:
            avg_depth = int(round(depth / num_lines))

    print(avg_depth)

if __name__ == "__main__":
    argv = sys.argv[1:]
    get_read_depths(argv[0])