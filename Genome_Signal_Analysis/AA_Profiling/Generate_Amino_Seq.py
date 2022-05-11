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

def add_intron(end, start, n, N, sequen, dstart, dend):
    if n >= N-1:
        if len(dstart) > 0:
            down_seq = sequen[int(dstart[0]):]
            down_seq = down_seq.replace('D', 'd')
            down_seq = down_seq.replace('I', 'i')
            return down_seq

        else: return 'f'*len(sequen[int(end[n]):])
    else:
        return 'i'*(int(start[n+1])- int(end[n]))



def formatting_file(n):


    aa_seq_file = open(AA_FILE_PATH+"\/region_group_" + str(n) + ".fasta", "w")

    print(n)
    for (seq_record_1, seq_record_2) in zip(
            SeqIO.parse(REGION_FILE_PATH + "\/region_group_" + str(n) + ".fasta", "fasta"),
            SeqIO.parse(SEQ_FILE_PATH + "\group_" + str(n) + ".fasta", "fasta")):
        region_header = seq_record_1.id
        seq_header = seq_record_2.id
        region = seq_record_1.seq[1000:]
        seq = seq_record_2.seq[1000:]

        region_array = np.array(list(region))
        region_array = np.where(region_array == 'E', 0, 1)
        new_region_arr = np.array(region_array, dtype=bool)
        masked_array = np.ma.masked_array(np.array(list(seq)), mask=new_region_arr)

        mrna = ''.join(np.ma.compressed(masked_array))
        aa_seq = Seq(mrna).translate(stop_symbol="x")

        # print(region_header)

        # print(seq_header)
        # print(len(aa_seq))

        region_header_ =region_header.split('|')
        if len(region_header_[1])>0:
            uutr_start = [int(s) for s in region_header_[1].split(';')]
        else:uutr_start = []

        if len(region_header_[2])>0:
            uutr_end = [int(s) for s in region_header_[2].split(';')]
        else:uutr_end=[]

        if len(region_header_[3])>0:
            dutr_start = [int(s) for s in region_header_[3].split(';')]
        else:dutr_start=[]

        if len(region_header_[4])>0:
            dutr_end = [int(s) for s in region_header_[4].split(';')]
        else: dutr_end=[]

        if len(region_header_[5])>0:
            exon_start = [int(s) for s in region_header_[5].split(';')]
        else:exon_start=[]

        if len(region_header_[6])>0:
            exon_end = [int(s) for s in region_header_[6].split(';')]
        else:exon_end = []

        first_index = 0
        last_index = 0
        aa_region = ''
        mod = 0

        for i in range(len(exon_start)):
            if mod != 0:
                exon_len = int(exon_end[i]) - int(exon_start[i])-(3-mod)
            else:
                exon_len = int(exon_end[i]) - int(exon_start[i])
            q, mod = divmod(exon_len, 3)
            # print(last_index)
            last_index = last_index + q
            # print('Quoteint and Remainder:', q, mod, ' Index Range: ',first_index, last_index)
            aa_region = aa_region + duplicate_aa(aa_seq[first_index:last_index], 3)

            if mod != 0:
                try:
                    aa_region = aa_region + duplicate_aa(aa_seq[last_index], mod)
                    aa_region = aa_region + add_intron(exon_end, exon_start, i, len(exon_start), seq_record_1.seq,
                                                       dutr_start, dutr_end)
                    aa_region = aa_region + duplicate_aa(aa_seq[last_index], 3 - mod)
                    last_index = last_index + 1
                    first_index = last_index

                except IndexError as error:
                    print(error, 'file_number', n)
                    print('erfd-----------------------', region_header_)
                    aa_region = aa_region + add_intron(exon_end, exon_start,i, len(exon_start),  seq_record_1.seq, dutr_start, dutr_end)


            else:
                first_index = last_index
                aa_region = aa_region + add_intron(exon_end, exon_start,i, len(exon_start),  seq_record_1.seq, dutr_start, dutr_end)

        if len(uutr_start)>0:
            up_seq = seq_record_1.seq[:uutr_end[-1]]
            up_seq = up_seq.replace('F', 'z')
            up_seq = up_seq.replace('U', 'u')
            up_seq = up_seq.replace('I', 'i')
            aa_region = up_seq + aa_region
        else:
            aa_region = 'z'*1000 + aa_region

        # print(aa_region)
        aa_seq_file.write(">" + region_header + "\n")
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