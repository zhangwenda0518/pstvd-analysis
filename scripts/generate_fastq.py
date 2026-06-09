import os
import sys
import re
import numpy as np
import random
import math


def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N',
                  'R': 'Y', 'M': 'K', 'W': 'W', 'S': 'S', 'Y': 'R', 'K': 'M',
                  'H': 'D', 'D': 'H', 'B': 'V', 'V': 'B'}
    return ''.join([complement[base] for base in dna[::-1]])


def parse_fasta(input_path):
    
    seq_id = ''
    seq_seq = ''
    
    with open(input_path, 'r') as infh:
        for buf in infh:
            buf = buf.replace('\n', '')
            
            if buf[0:1] == '>':
                if seq_id != '' and seq_seq != '':
                    seq_id = ''
                    seq_seq = ''
                
                seq_id = buf[1:].split(' ')[0]
            
            else:
                seq_seq = seq_seq + buf
    
    return [seq_id, seq_seq]




def simulate_shortseq(seq, prefix, seq_len, x_coverage, seq_len_type, init_seed):
    random.seed(init_seed)
    
    read_names = []
    read_seqs  = []
    
    n_seq = int(len(seq) * x_coverage / seq_len)
    
    for i in range(n_seq):
        read_len = 0
        if seq_len_type == 'FIXED':
            read_len = seq_len
        elif seq_len_type == 'POISSON':
            j = 0
            while read_len < 18 or 26 < read_len:
                read_len = np.random.poisson(lam=seq_len)
                j = j + 3
        else:
            raise ValueError('use `FIXED` or `POISSON`.')
        
        read_start = np.random.randint(0, len(seq))
        
        ## slice the short seqeuneces from circle RNA
        seq_extend = seq + seq
        read_seq = seq_extend[read_start:(read_start + read_len)]
        
        strand = 'P'
        
        if random.random() > 0.5:
            strand = 'N'
            read_seq = reverse_complement(read_seq)
        
        read_names.append('{0}.{1:0=6}_{2}'.format(prefix, i, strand))
        read_seqs.append(read_seq)
    
    return (read_names, read_seqs)
    
    

def generate_fastq(input_path, output_path, seq_len, x_coverage, seq_len_type):
    
    _id, _seq = parse_fasta(input_path)
        
    with open(output_path, 'w') as outfh:
        read_names, read_seqs = simulate_shortseq(_seq, _id, seq_len, x_coverage, seq_len_type, init_seed = 202020 + hash(_id) % 100000)
        
        for read_name, read_seq in zip(read_names, read_seqs):
            outfh.write('@{0}\n{1}\n+\n{2}\n'.format(read_name, read_seq, 'I' * len(read_seq)))



if __name__ == '__main__':
    
    project_path = os.path.dirname(os.path.abspath(__file__))
    
    if len(sys.argv) != 5:
        raise ValueError('Usage: python $this input.fa output.fq 16 FIXED')
    
    
    input_fasta = sys.argv[1]
    output_fastq = sys.argv[2]
    seq_len = int(sys.argv[3])
    seq_len_type = sys.argv[4]
    
    x_coverage = int(os.environ.get("PSTVD_COVERAGE", 10000))
    
    generate_fastq(input_fasta, output_fastq, seq_len, x_coverage, seq_len_type)
    
    


