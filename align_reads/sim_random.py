import numpy as np

rndg = np.random.default_rng(454)

# base bias
#tcount = 320
#gcount = 929
#ccount = 645
#acount = 2054

tcount = 100
gcount = 100
ccount = 100
acount = 100

total_count = acount + ccount + gcount + tcount

tcut = 1.0 - (tcount)/total_count
gcut = 1.0 - (tcount+gcount)/total_count
ccut = 1.0 - (tcount+gcount+ccount)/total_count

def gen_read(readlen):
    read = ''
    for i in range(readlen):
        rnum = rndg.random()
        if rnum > tcut:
            b = 'T'
        elif rnum > gcut:
            b = 'G'
        elif rnum > ccut:
            b = 'C'
        else:
            b = 'A'
        read += b
    return read

outfile = open('reads.random.txt', 'w')
num_reads = 3833
readlen = 13

for i in range(num_reads):
    read = gen_read(readlen)
    outfile.write('%s\n' % read)
outfile.close()

