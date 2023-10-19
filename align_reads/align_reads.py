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

def bases2vals(bases):
    valdict = {'A':'0', 'C':'1', 'G':'2', 'T':'3'}
    vals = ''
    for base in bases:
        vals += valdict[base]
    return vals

# generates an array of size 4^hash_len and each entry contains a list of genome positions matching the hash
def generate_hashlist(genome, hash_len):
    hash_list = []
    genome_len = len(genome)
    num_hashs = int(math.pow(4,hash_len))

    for i in range(num_hashs):
        hash_list.append([])

    genome_position = 0
    while genome_position < (genome_len-hash_len):
        index = int(bases2vals(genome[genome_position:genome_position+hash_len]), 4)
        hash_list[index].append(genome_position)
        genome_position += 1

    return hash_list

def map_read(genome, read, hash_list, hash_len):
    # for each genome position in the hash list, compute edit distance for read, return the position & edit distance of the best match
    hash_index = int(bases2vals(read[:hash_len]), 4)
    bestdist = -1
    bestpos = -1
    for pos in hash_list[hash_index]:
        readlen = len(read)
        dist = editdistance.eval(genome[pos:pos+readlen], read)
        if bestdist == -1 or dist < bestdist:
            bestdist = dist
            bestpos = pos

    return bestpos,bestdist

# set some defaults
ref_name = 'phix174.fasta'
verbose = 0
in_filename = None
out_filename = None
hash_len = 4

# parce cmd-line args
argcc = 1
argc = len(sys.argv)
while argcc < argc:
    if sys.argv[argcc] == '--ref':
        argcc += 1
        ref_name = sys.argv[argcc]
    if sys.argv[argcc] == '--in':
        argcc += 1
        in_filename = sys.argv[argcc]
    if sys.argv[argcc] == '--out':
        argcc += 1
        out_filename = sys.argv[argcc]
    if sys.argv[argcc] == '--hl':
        argcc += 1
        hash_len = int(sys.argv[argcc])
    if sys.argv[argcc] == '-v':
        verbose += 1
    argcc += 1

# load up our reads
if in_filename is None:
    print('missing --in readfile?')
    exit(0)
reads = []
with open(in_filename) as f:
    lines = f.readlines()
    for line in lines:
        reads.append(line.rstrip())

# load up the fasta file
ref,description = load_fasta(ref_name)
if verbose > 0:
    print('loaded ref: %s\nlength: %d\n' % (description, len(ref)))

if verbose > 1:
    print('first few reads:')
    for r in range(20):
        print('%s' % reads[r])

# generate hash list from genome and map reads
print('generating hash list...')
hash_list = generate_hashlist(ref, hash_len)

outfile = None
if out_filename:
    outfile = open(out_filename, 'w') 

print('mapping reads...')
for read in reads:
    pos,dist = map_read(ref, read, hash_list, hash_len)
    if pos > -1:
        outtxt = 'read: %s  ref: %s  pos: %d  dist: %d' % (read, ref[pos:pos+len(read)], pos, dist)
        if outfile:
            outfile.write('%s\n' % outtxt)
        else:
            print(outtxt)

if outfile:
    outfile.close()

