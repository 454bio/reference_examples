import numpy as np
import math
import sys

def qualscore(readlen, num_errors):
    qual_score = -10.0 * math.log10(num_errors/readlen)
    return qual_score

def scoremin(seq, ref, minlen):
    bestq = 0
    best_q_len = 0
    perfect_len = 0

    errors = 0
    ls = len(seq)
    lr = len(ref)
    if lr < ls:
        ls = lr

    q = 0

    for i in range(minlen):
        if seq[i] != ref[i]:
            errors += 1
        if perfect_len == 0 and errors == 1:
            perfect_len = i
    if errors > 0:
        q = qualscore(minlen, errors)
        bestq = q
        best_q_len = minlen
    else:
        perfect_len = minlen

    for i in range(minlen,ls):
        if seq[i] != ref[i]:
            errors += 1
        if errors > 0:
            q = qualscore(i+1, errors)
        else:
            q = 0
        if perfect_len == 0 and errors == 1:
            perfect_len = i
        if q > bestq:
            bestq = q
            best_q_len = i+1

    return (bestq,best_q_len,perfect_len)

def top_n_by_len(N, L):
    forced_info = []
    for i in range(len(seqlist)):
        # start calculating the Q-Score at length L (our minimum) and continue to the end of the read to get the best Q-Score
        forced_info.append(scoremin(seqlist[i], reflist[i], L))
    forced_info = np.array(forced_info)

    top_by_q_i = np.flip(np.argsort(forced_info[:,0]))
    top_by_q = forced_info[top_by_q_i]
    print('%d %dQ%.2f' % (N, L, np.mean(top_by_q[:N,0])))

# defaults
filename = '20231009_S0495.out'

# process cmd-line args
argc = len(sys.argv)
argcc = 1
while argcc < argc:
    if sys.argv[argcc] == '-f':
        argcc += 1
        filename = sys.argv[argcc]
    argcc += 1

# load up the results file
seqlist = []
reflist = []
with open(filename) as f:
    lines = f.readlines()
    for line in lines:
        t = line.split(' ')
        try:
            seq = t[7]
            ref = t[12]
            seqlist.append(seq)
            reflist.append(ref)
        except:
            pass

# calculate stats based on min read length and top N for various combinations
top_n_by_len(100, 6)
top_n_by_len(50, 10)
top_n_by_len(100, 10)
top_n_by_len(100, 12)

