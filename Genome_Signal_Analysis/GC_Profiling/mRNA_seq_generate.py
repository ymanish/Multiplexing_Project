'''The script will read the unspliced transcripts and return the mrna sequence and region information in a single file'''

from Bio import SeqIO
import concurrent.futures
import time
from Initialise_GC_profiling import *

def generate_mrna_region(header, region_sequ, gene_sequ):

    ID = header.split('|')[0]

    if len(header.split('|')[1]) >0:
        UUTR_loc_start = header.split('|')[1].split(';')
        UUTR_loc_stop = header.split('|')[2].split(';')
    else:
        UUTR_loc_start = []
        UUTR_loc_stop = []

    if len(header.split('|')[3])>0:
        DUTR_loc_start = header.split('|')[3].split(';')
        DUTR_loc_stop = header.split('|')[4].split(';')
    else:
        DUTR_loc_start = []
        DUTR_loc_stop = []

    if len(header.split('|')[5])>0:
        exon_loc_start = header.split('|')[5].split(';')
        exon_loc_stop = header.split('|')[6].split(';')
    else:
        exon_loc_start = []
        exon_loc_stop = []

    new_region = region_sequ[:1000]
    mrna_seq = gene_sequ[:1000]
    # new_region = ''
    # mrna_seq = ''

    UUTR_len = 1000
    if len(UUTR_loc_start) > 0:

        for s, e in zip(UUTR_loc_start, UUTR_loc_stop):
            new_region = new_region + region_sequ[int(s):int(e)]
            mrna_seq = mrna_seq + gene_sequ[int(s):int(e)]
            UUTR_len = UUTR_len + (int(e)-int(s))


    exon_len = UUTR_len
    if len(exon_loc_start) > 0:

        for s, e in zip(exon_loc_start, exon_loc_stop):
            new_region = new_region + region_sequ[int(s):int(e)]
            mrna_seq = mrna_seq + gene_sequ[int(s):int(e)]
            # print(int(e) - int(s))
            exon_len = exon_len + (int(e) - int(s))


    DUTR_len = exon_len

    if len(DUTR_loc_start) > 0:
        for s, e in zip(DUTR_loc_start, DUTR_loc_stop):
            new_region = new_region + region_sequ[int(s):int(e)]
            mrna_seq = mrna_seq + gene_sequ[int(s):int(e)]
            # print(int(e) - int(s))
            DUTR_len = DUTR_len + (int(e) - int(s))


    mrna_seq = mrna_seq + gene_sequ[len(gene_sequ)-1000:]
    new_region = new_region + region_sequ[len(region_sequ)-1000:]


    new_header = str(UUTR_len) + '|'+ str(exon_len)+ '|'+str(DUTR_len)
    return new_region, mrna_seq, ID, new_header

def main(n):

    mrna_region_file = open(REGION_FILE_PATH+"\/region_group_" + str(n) + ".fasta", "w")
    mrna_seq_file = open(SEQ_FILE_PATH+"\/group_" + str(n) + ".fasta", "w")

    for (seq_record_1, seq_record_2) in zip(
            SeqIO.parse(pre_mrna_REGION_FILE+"\/region_group_"+str(n)+".fasta", "fasta"),
            SeqIO.parse(pre_mrna_SEQ_FILE+"\group_"+str(n)+".fasta", "fasta")):

        region_header = seq_record_1.id
        region = seq_record_1.seq

        sequence = seq_record_2.seq
        new_region, new_seq, id, n_header = generate_mrna_region(region_header, region, sequence)
        print(id)

        mrna_region_file.write(">" + id + '|' + n_header + "\n")
        mrna_region_file.write(str(new_region) + "\n")

        mrna_seq_file.write(">" + id + '|' + n_header + "\n")
        mrna_seq_file.write(str(new_seq) + "\n")

    mrna_region_file.close()
    mrna_seq_file.close()

    return n

if __name__ == "__main__":

    start = time.perf_counter()

    with concurrent.futures.ProcessPoolExecutor() as executor:
        iter_seq = range(1, Total_Input_Files+1)
        pool = [executor.submit(main, n=i) for i in iter_seq]
        for i in concurrent.futures.as_completed(pool):
            print(f'Return Value: {i.result()}')

    end = time.perf_counter()
    print(f'Finished in {round(end - start, 2)} second(s)')
