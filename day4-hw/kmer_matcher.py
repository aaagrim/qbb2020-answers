#!/usr/bin/env python3

from fasta_iterator_class import FASTAReader


kmers = {}

our_seq = {}

k = 11

for seq_id, sequence in FASTAReader(open('droYak2_seq.fa')):    
    for i in range(0, len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        kmers.setdefault(kmer,[])
        kmers[kmer].append(i)

for seq_id, sequence in FASTAReader(open('subset.fa')):
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