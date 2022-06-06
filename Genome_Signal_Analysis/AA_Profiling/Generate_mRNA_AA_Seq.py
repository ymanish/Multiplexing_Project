import numpy as np
from Bio import SeqIO
import concurrent.futures
import time
from Initialise_AA_profiling import *
from Bio.Seq import Seq


def duplicate_aa(aa, m):
    aa_ = ''
    for a in aa:
        aa_ = aa_ + a*m
    return aa_

def formatting_file(n):


    aa_seq_file = open(mrna_AA_FILE_PATH+"\/region_group_" + str(n) + ".fasta", "w")

    print(n)
    for (seq_record_1, seq_record_2) in zip(
            SeqIO.parse(mrna_REGION_FILE_PATH + "\/region_group_" + str(n) + ".fasta", "fasta"),
            SeqIO.parse(mrna_SEQ_FILE_PATH + "\group_" + str(n) + ".fasta", "fasta")):

        region_header = seq_record_2.id.split('|')

        q, mod = divmod(int(region_header[2])-int(region_header[1]), 3)

        if mod!=0:
            # raise ValueError('The length of the CDS is not a multiple of three '+region_header[0])
            mrna = seq_record_2.seq[int(region_header[1]):int(region_header[2])]
            mrna = mrna + 'N'*(3-mod)
            aa_seq = Seq(mrna).translate(stop_symbol="x")
            aa_seq = duplicate_aa(aa_seq, 3)
            aa_seq = aa_seq[:-(3-mod)]
        else:
            mrna = seq_record_2.seq[int(region_header[1]):int(region_header[2])]
            aa_seq = Seq(mrna).translate(stop_symbol="x")
            aa_seq = duplicate_aa(aa_seq, 3)


        aa_region = 'z'*1000
        aa_region = aa_region+'u'*(int(region_header[1])-1000)
        aa_region = aa_region + aa_seq
        aa_region = aa_region +'d'*(int(region_header[3])-int(region_header[2]))
        aa_region = aa_region +'f'*(1000)

        assert len(aa_region) == len(seq_record_2.seq), "length of the sequence and region does not match for " + region_header[0]

        aa_seq_file.write(">" + '|'.join(region_header) + "\n")
        aa_seq_file.write(str(aa_region) + "\n")
    aa_seq_file.close()

    return None

if __name__ == "__main__":

    start = time.perf_counter()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
        pool = [executor.submit(formatting_file, n=i) for i in iter_seq]
        for i in concurrent.futures.as_completed(pool):
            print(f'Return Value: {i.result()}')

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')