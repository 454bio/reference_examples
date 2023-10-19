import sys
import numpy as np
import math
import editdistance

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

# generates normally distributed readlengths of reads covering a uniform genome at the specified depth
def generate_reads(genome, mean_length, coverage_mult):
    reads = []
    genome_len = len(genome)
    num_reads = int(coverage_mult*genome_len/mean_length)
    rnd = np.random.default_rng(454)
    for r in range(num_reads):
        # pick a uniform random start location in the genome
        start = rnd.integers(genome_len-mean_length)
        length = mean_length + int(mean_length*0.2*rnd.normal(0, 1))
        reads.append(genome[start:start+length])
    return reads

def generate_errors(read):
    bases = 'ACGT'
    pos = 4 # start at the 5th position assuming the first 4 bases are always 100% correct
    readlen = len(read)
    rnd = np.random.default_rng(1134)
    read = list(read)
    while pos < readlen:
        if rnd.random() > 0.8: # 20% chance of error
            newbase = rnd.integers(4) # range is 0 to 3 inclusive
            read[pos] = bases[newbase]
        pos += 1
    return ''.join(read)

# set some defaults
ref_name = 'phix174.fasta'
verbose = 0
out_filename = None
num_cycles = 30

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
    if sys.argv[argcc] == '-v':
        verbose += 1
    argcc += 1

# load up the fasta file
ref,description = load_fasta(ref_name)
if verbose > 0:
    print('loaded ref: %s\nlength: %d\n' % (description, len(ref)))

# generate some random reads from it
reads = generate_reads(ref, 200, 3) # mean 200bp at 3x coverage

# truncate the reads to our num cycles
for i in range(len(reads)):
    reads[i] = reads[i][:num_cycles]

# add some random errors
for i in range(len(reads)):
    reads[i] = generate_errors(reads[i])

if verbose > 1:
    print('first few reads:')
    for r in range(20):
        print('%s' % reads[r])

# dump the reads to a file
if out_filename is not None:
    with open(out_filename, 'w') as f:
        for read in reads:
            f.write(read)
            f.write('\n')

