from itertools import islice
import sys
import re
ENTRY_LEN = 4
CIGAR_LOC = 5
READS_LOC = 9
QUALS_LOC = 10
READID_LOC = 0
NEXT_ID_LOC = 4
ALPHABET = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']

def split_bam(output, f_id):
    new_fastq = open(output + "/" + f_id + "_splitgap.fastq", "w")
    #unsplitables = open(f_id + "_no_split.fastq", "w")
    oldfq_id = 0
    with open(output + "/" + f_id + "_bam.txt") as bam:
        with open(output + "/" + f_id + "_old.fastq") as old_fastq:
            #split old fqs into array format per line
            oldfq = old_fastq.read().split()
            for line in bam:
                arr = str(line).split()
                read_id, cigar, seq, quals = parse_line(arr) 
                fq_read_id = oldfq[oldfq_id]
                new_seq, new_quals = split_reads_and_quals(cigar, seq, quals)
                #paired end in the fastq are swapped sometimes in the bam file, so as long as read_id's are the same print whatever paired end sequence comes up first
                if read_id in fq_read_id:
                    write_to_fastq(new_fastq, fq_read_id, new_seq, new_quals, oldfq_id)
                    oldfq_id += NEXT_ID_LOC
                #read was not in the fastq file
                else:
                    continue
    new_fastq.close()

def write_to_fastq(fastq, read_id, new_seq, new_quals, fq_id):
    for i in range(0, len(new_seq)):
        fastq.write(get_next_id(read_id, i) + "\n")
        fastq.write(new_seq[i] + "\n")
        fastq.write("+\n")
        fastq.write(new_quals[i] + "\n")

def get_next_id(read_id, next_let):
    return read_id + ALPHABET[next_let]

def parse_line(arr):
    return arr[READID_LOC], arr[CIGAR_LOC], arr[READS_LOC], arr[QUALS_LOC]

def split_reads_and_quals(cigar, reads, quals):
    if cigar == "*":
        return [reads], [quals]
    new_reads = []
    new_quals = []
    #regex matches any amount of numbers followed by one capital letter
    matches = re.findall(r'[0-9]+[A-Z]{1}', cigar)
    curr_num = 0
    temp_idx = 0
    let = ""
    for match in matches:
        num = int(match[:-1])
        let = match[-1]
        #append next reads and quals if not an M or an I
        if let != 'M' or let != 'I' and curr_num != 0:
            add_next_seq(temp_idx, curr_num, new_reads, reads, new_quals, quals)
            temp_idx += curr_num
            curr_num = 0
        else:
            curr_num += num
    #append last reads/quals if M or I
    if let == 'M' or let == 'I':
        add_next_seq(temp_idx, curr_num, new_reads, reads, new_quals, quals)   
        temp_idx += curr_num 

    return new_reads, new_quals

def add_next_seq(temp_idx, curr_num, new_reads, reads, new_quals, quals):
    next_idx = temp_idx+curr_num
    new_reads.append(reads[temp_idx:next_idx])
    new_quals.append(quals[temp_idx:next_idx])

if __name__ == "__main__":
    argv = sys.argv[1:]
    split_bam(argv[0], argv[1])
