import sys
import numpy as np


filename = 'stats.out'

argcc = 1
argc = len(sys.argv)
while argcc < argc:
    if sys.argv[argcc] == '-i':
        argcc += 1
        filename = sys.argv[argcc]
    argcc += 1

qualScores = []
perfectScores = []
ie = []
cf = []
dr = []
perfectHist = [0]*13

with open(filename) as f:
    lines = f.readlines()
    for line in lines:
        tokens = line.strip().split(' ')
        if len(tokens) >= 19:
            qualScores.append(float(tokens[16]))

            dist = int(tokens[18])
            perfectHist[dist] += 1
            if dist >= 6:
                perfectScores.append(float(tokens[16]))

            err = float(tokens[5])
            if err < 3.0:
                phase = tokens[3].split('/')
                ie.append(float(phase[0]))
                cf.append(float(phase[1]))
                dr.append(float(phase[2]))


qualScores = np.sort(np.array(qualScores))
print('avg qual score for top 50: %.3f' % np.mean(qualScores[-50:]))
print('avg qual score for top 100: %.3f' % np.mean(qualScores[-100:]))

perfectScores = np.array(perfectScores)
print('%d 6 Q %.3f' % (len(perfectScores), np.mean(perfectScores)))

ie = np.array(ie)
cf = np.array(cf)
dr = np.array(dr)

print('%d samples have average ie/cf/dr: %.3f/%.3f/%.3f' % (len(ie), np.mean(ie), np.mean(cf), np.mean(dr)))
print('%d samples have median ie/cf/dr: %.3f/%.3f/%.3f' % (len(ie), np.median(ie), np.median(cf), np.median(dr)))

print('perfect length distribution: %s' % perfectHist)

