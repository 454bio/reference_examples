basecall.py
a reference example from the jmr twitter dataset
updated to include a 4-dye variable rate model

Usage:

basic run:
python3 ./basecall.py -ie 0.09 -cf 0.065 -dr 0.02 -plots
output:
cycles: 40
basecalls: CTGACTGACTAACGTCAGTCCATGCACTCACGTTTAACTG
cumulative error: 0.704400
also creates a couple of plots to show original input signals per cycle, predicted signals, error, signal loss

grid-search
python3 ./basecall.py -ie 0.09 -cf 0.065 -dr 0.02 -plots -grid
outputs base calls for each cf/ie/dr combination, tracks total error per combination, and selects the lowest error rate

4-dye run:
python3 ./basecall.py -ie 0.09 -cf 0.065 -dr 0.02 -plots -model 4dye
note that the 4 dye version is hard-coding a 'G' at 1.5 higher cf rate as an example
the rates should all be relative to some baseline, and the grid search is not designed to work with this version

