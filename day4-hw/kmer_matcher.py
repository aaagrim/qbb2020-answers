#!/usr/bin/env python3

import sys

from fasta_iterator_class import FASTAReader

target_file = sys.argv[1]
query_file = sys.argv[2]

k = int(sys.argv[3])

kmers = {}

our_seq = {}


for seq_id, sequence in FASTAReader(open(query_file)):    
    for i in range(0, len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        kmers.setdefault(kmer,[])
        kmers[kmer].append(i)

for seq_id, sequence in FASTAReader(open(target_file)):
    for i in range(0, len(sequence) - k + 1):
        t_kmer = sequence[i:i + k]
        our_seq.setdefault(t_kmer,{})
        our_seq[t_kmer].setdefault(seq_id,[])
        our_seq[t_kmer][seq_id].append(i) 
    
for q_seq in kmers.keys():
    if q_seq in our_seq.keys():
        for nucs in our_seq[q_seq]:
            print(nucs,our_seq[q_seq][nucs], kmers[q_seq], q_seq)
            
if __name__ == "__main__":
    reader = FASTAReader(sys.stdin)

    for seq_id, sequence in reader:
        print(seq_id, sequence)