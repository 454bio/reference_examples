import sys
import numpy as np
import math
import editdistance

# global random number generator
rndg = np.random.default_rng(1134)

# simple function to load a fasta file
def load_fasta(filename):
    genome = ''
    description = ''
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.rstrip()
            if line[0] == '>':
                description = line[1:]
            else:
                genome += line

    return genome, description

def reverse_comp(read):
    rcomp = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    rread = ''
    for c in read[::-1]:
        rread += rcomp[c]
    return rread

# generates normally distributed readlengths of reads covering a uniform genome at the specified depth
def generate_reads(genome, mean_length, num_reads):
    info = []
    genome_len = len(genome)
    rnd = np.random.default_rng(454)
    for r in range(num_reads):
        # pick a uniform random start location in the genome
        start = rnd.integers(genome_len-mean_length)
        length = mean_length + int(mean_length*0.2*rnd.normal(0, 1))
        read = genome[start:start+length]
        # determine reverse complement or normal 50/50 chance
        rcomp = False
        if rnd.normal(0, 1) > 0.0:
            read = reverse_comp(read)
            rcomp = True
        info.append({'read':read, 'rcomp':rcomp, 'start':start, 'length':length})
    return info

def generate_errors(read, err_rate):
    bases = 'ACGT'
    readlen = len(read) - 4 # we will assume the first 4 bases are always 100% correct
    pos = 4
    read = list(read)
    while pos < readlen:
        if rndg.random() > (1.0-err_rate): # % chance of error
            newbase = rndg.integers(4) # range is 0 to 3 inclusive
            read[pos] = bases[newbase]
        pos += 1
    return ''.join(read)

# set some defaults
ref_name = 'phix174.fasta'
verbose = 0
out_filename = None
num_cycles = 30
err_rate = 0.2
num_reads = 10000

# parce cmd-line args
argcc = 1
argc = len(sys.argv)
while argcc < argc:
    if sys.argv[argcc] == '--ref':
        argcc += 1
        ref_name = sys.argv[argcc]
    if sys.argv[argcc] == '--out':
        argcc += 1
        out_filename = sys.argv[argcc]
    if sys.argv[argcc] == '--cycles':
        argcc += 1
        num_cycles = int(sys.argv[argcc])
    if sys.argv[argcc] == '--err_rate':
        argcc += 1
        err_rate = float(sys.argv[argcc])
    if sys.argv[argcc] == '--num_reads':
        argcc += 1
        num_reads = int(sys.argv[argcc])
    if sys.argv[argcc] == '-v':
        verbose += 1
    argcc += 1

# load up the fasta file
ref,description = load_fasta(ref_name)
if verbose > 0:
    print('loaded ref: %s\nlength: %d\n' % (description, len(ref)))

# generate some random reads from it
info = generate_reads(ref, 200, num_reads) # mean 200bp at 3x coverage

# truncate the reads to our num cycles
reads = []
for i in range(len(info)):
    reads.append(info[i]['read'][:num_cycles])

# add some random errors
for i in range(len(reads)):
    reads[i] = generate_errors(reads[i], err_rate)

if verbose > 1:
    print('raw read info:')
    for r in range(len(reads)):
        print('raw read: %s read: %s start: %d len: %d rcomp: %s' % (info[r]['read'][:num_cycles], reads[r], info[r]['start'], info[r]['length'], info[r]['rcomp']))

# dump the reads to a file
if out_filename is not None:
    with open(out_filename, 'w') as f:
        for read in reads:
            f.write(read)
            f.write('\n')

