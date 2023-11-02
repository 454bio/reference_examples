import sys

in_filename = 'reads.out'
out_filename = 'reads.txt'

argcc = 1
argc = len(sys.argv)
while argcc < argc:
    if sys.argv[argcc] == '-i' or sys.argv[argcc] == '-f':
        argcc += 1
        in_filename = sys.argv[argcc]
    if sys.argv[argcc] == '-o':
        argcc += 1
        out_filename = sys.argv[argcc]
    argcc += 1

out_file = open(out_filename, 'w')

with open(in_filename) as f:
    lines = f.readlines()
    for line in lines:
        tokens = line.strip().split(' ')
        if len(tokens) >= 8:
            out_file.write('%s\n' % tokens[7])

out_file.close()

