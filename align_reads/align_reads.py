import sys
import numpy as np
import math
import editdistance
import sam_utils

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
    bestdist = -1
    bestpos = -1
    readlen = len(read)
    hash_index = int(bases2vals(read[:hash_len]), 4) # cool python way to convert a base-4 string to a base-10 integer
    for pos in hash_list[hash_index]:
        dist = editdistance.eval(genome[pos:pos+readlen], read)
        if bestdist == -1 or dist < bestdist:
            bestdist = dist
            bestpos = pos

    return bestpos,bestdist

def eval_read(read, ref_read, num_errors):
    err = 0
    readlen = len(read)
    maxlen = 0
    while maxlen < readlen:
        if read[maxlen] != ref_read[maxlen]:
            err += 1
        if err > num_errors:
            break
        maxlen += 1
    return maxlen

# set some defaults
ref_name = 'phix174.fasta'
verbose = 0
in_filename = None
out_filename = None
hash_len = 4

# parse cmd-line args
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
ref_r = reverse_comp(ref)
hash_list_r = generate_hashlist(ref_r, hash_len)

outfile = None
if out_filename:
    outfile = open(out_filename, 'w') 
    sam_filename = out_filename.split('.')[0] + '.sam'
else:
    sam_filename = 'out.sam'

# some stats we will track while mapping
num12Q10 = 0

# map reads and write sam file
sam = sam_utils.SamUtils(sam_filename)

print('mapping reads...')
for i, read in enumerate(reads):
    rcomp = False
    # we map to both the forward and the reverse complement of the read, to see which is the best match
    pos,dist = map_read(ref, read, hash_list, hash_len)
    pos_r,dist_r = map_read(ref_r, read, hash_list_r, hash_len)
    if dist_r < dist:
        dist = dist_r
        pos = pos_r
        rcomp = True
    if pos > -1:
        ref_read = ref_r[pos:pos+len(read)] if rcomp else ref[pos:pos+len(read)]
        if len(ref_read) == len(read):
            outtxt = 'read: %s  ref: %s  rcomp: %s  pos: %d  dist: %d' % (read, ref_read, rcomp, pos, dist)
            if outfile:
                outfile.write('%s\n' % outtxt)
            else:
                print(outtxt)
            read_name = 'read_' + str(i)
            sam.AddRead(read, pos, ref_read, read_name)

            # calc some stats on the read
            len1Error = eval_read(read, ref_read, 1) # the max length with only 1 error
            if len1Error >= 12:
                num12Q10 += 1
        else:
            print('failed to align read: %s to valid reference position' % read)

if outfile:
    outfile.close()

print('%d 12Q10' % num12Q10)

