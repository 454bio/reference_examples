python3 sim_reads.py --ref phix174.fasta --out reads.txt --cycles 25
python3 align_reads.py --ref phix174.fasta --in reads.txt --out aligned_reads.txt

