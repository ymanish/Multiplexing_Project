from Bio import SeqIO
import copy
import concurrent.futures
import time
from Initialise_GC_profiling import *

#Eukaryote Sequence fasta header info and order:
#"ensembl_gene_id"
#"ensembl_gene_id_version"
#"ensembl_transcript_id"
#"ensembl_transcript_id_version"
#"5_utr_start"
#"5_utr_end"
#"3_utr_start"
#"3_utr_end"
#"exon_chrom_start"
#"exon_chrom_end"
#"transcription_start_site"
#"strand"
#"rank"
#"chromosome_name"


def genomic_positions_eukaryote(strand, header):
    if strand == 1:
        flag = False
    else:
        flag = True

    GENE_ID = header[0]
    TRANSCRIPT_ID = header[2]

    if len(header[4]) > 0:
        FIVE_UTR_START = [int(i) for i in header[4].split(';')]
        FIVE_UTR_START = sorted(FIVE_UTR_START, reverse=flag)
    else:
        FIVE_UTR_START = []
    if len(header[5]) > 0:
        FIVE_UTR_END = [int(i) for i in header[5].split(';')]
        FIVE_UTR_END = sorted(FIVE_UTR_END, reverse=flag)
    else:
        FIVE_UTR_END = []

    if len(header[6]) > 0:
        THREE_UTR_START = [int(i) for i in header[6].split(';')]
        THREE_UTR_START = sorted(THREE_UTR_START, reverse=flag)
    else:
        THREE_UTR_START = []

    if len(header[7]) > 0:
        THREE_UTR_END = [int(i) for i in header[7].split(';')]
        THREE_UTR_END = sorted(THREE_UTR_END, reverse=flag)
    else:
        THREE_UTR_END = []

    EXON_STARTS = [int(i) for i in header[8].split(';')]  # array
    EXON_ENDS = [int(i) for i in header[9].split(';')]  # array

    TSS = int(header[10])
    # STRAND=int(header[11])
    EXON_RANK = [int(i) for i in header[12].split(';')]  # array
    CHR = header[13]

    return GENE_ID, TRANSCRIPT_ID, FIVE_UTR_START, FIVE_UTR_END, THREE_UTR_START, THREE_UTR_END, EXON_STARTS, EXON_ENDS, TSS, EXON_RANK, CHR

def genomic_positions_plant(strand, header):
    if strand == 1:
        flag = False
    else:
        flag = True

    GENE_ID = header[0]
    TRANSCRIPT_ID = header[1]

    if len(header[2]) > 0:
        FIVE_UTR_START = [int(i) for i in header[2].split(';')]
        FIVE_UTR_START = sorted(FIVE_UTR_START, reverse=flag)
    else:
        FIVE_UTR_START = []
    if len(header[3]) > 0:
        FIVE_UTR_END = [int(i) for i in header[3].split(';')]
        FIVE_UTR_END = sorted(FIVE_UTR_END, reverse=flag)
    else:
        FIVE_UTR_END = []

    if len(header[4]) > 0:
        THREE_UTR_START = [int(i) for i in header[4].split(';')]
        THREE_UTR_START = sorted(THREE_UTR_START, reverse=flag)
    else:
        THREE_UTR_START = []

    if len(header[5]) > 0:
        THREE_UTR_END = [int(i) for i in header[5].split(';')]
        THREE_UTR_END = sorted(THREE_UTR_END, reverse=flag)
    else:
        THREE_UTR_END = []

    EXON_STARTS = [int(i) for i in header[6].split(';')]  # array
    EXON_ENDS = [int(i) for i in header[7].split(';')]  # array

    TSS = int(header[8])
    # STRAND=int(header[9])

    EXON_RANK = [int(i) for i in header[10].split(';')]  # array
    CHR = header[11]

    return GENE_ID, TRANSCRIPT_ID, FIVE_UTR_START, FIVE_UTR_END, THREE_UTR_START, THREE_UTR_END, EXON_STARTS, EXON_ENDS, TSS, EXON_RANK, CHR


def positive_strand_positions(FIVE_UTR_START, FIVE_UTR_END, THREE_UTR_START, THREE_UTR_END, TSS):
    uutr_starts = []
    uutr_stops = []

    dutr_starts = []
    dutr_stops = []

    for s, e in zip(FIVE_UTR_START, FIVE_UTR_END):
        uutr_starts.append(1000 + s - TSS)
        uutr_stops.append(1000 + (e - TSS) + 1)


    for s, e in zip(THREE_UTR_START, THREE_UTR_END):
        dutr_starts.append(1000 + s - TSS)
        dutr_stops.append(1000 + (e - TSS) + 1)

    return uutr_starts, uutr_stops, dutr_starts, dutr_stops


def negative_strand_positions(FIVE_UTR_START, FIVE_UTR_END, THREE_UTR_START, THREE_UTR_END, TSS):
    uutr_starts = []
    uutr_stops = []

    dutr_starts = []
    dutr_stops = []

    for s, e in zip(FIVE_UTR_START, FIVE_UTR_END):
        uutr_starts.append(1000 + TSS - e)
        uutr_stops.append(1000 + (TSS - s) + 1)

    for s, e in zip(THREE_UTR_START, THREE_UTR_END):
        dutr_starts.append(1000 + TSS - e)
        dutr_stops.append(1000 + (TSS - s) + 1)

    return uutr_starts, uutr_stops, dutr_starts, dutr_stops


def exon_rank(EXON_RANK, EXON_STARTS, EXON_ENDS, TSS, strand):
    EXON_STARTS = [i for _, i in sorted(zip(EXON_RANK, EXON_STARTS))]
    EXON_ENDS = [i for _, i in sorted(zip(EXON_RANK, EXON_ENDS))]

    #     G_CDS_START=[i for _,i in sorted(zip(CDS_START, G_CDS_START))]
    #     G_CDS_END=[i for _,i in sorted(zip(CDS_END, G_CDS_END))]

    exon_start = []
    exon_end = []

    if strand == 1:
        for s, e in zip(EXON_STARTS, EXON_ENDS):
            exon_start.append(1000 + s - TSS)
            exon_end.append(1000 + e - TSS + 1)

    else:

        for s, e in zip(EXON_STARTS, EXON_ENDS):
            exon_start.append(1000 + TSS - e)
            exon_end.append(1000 + TSS - s + 1)

    return exon_start, exon_end


def region_header(exon_start, exon_end, uutr_starts, uutr_stops, dutr_starts, dutr_stops, strand, chr):

    mexon_start = copy.deepcopy(exon_start)
    mexon_end = copy.deepcopy(exon_end)

    if len(uutr_starts) > 0 and len(dutr_starts) > 0:

        mexon_start[len(uutr_stops) - 1] = uutr_stops[-1]
        mexon_end[-len(dutr_starts)] = dutr_starts[0]

        mexon_start = mexon_start[len(uutr_stops) - 1:len(exon_start) - len(dutr_starts) + 1]
        mexon_end = mexon_end[len(uutr_stops) - 1:len(exon_start) - len(dutr_starts) + 1]

    elif len(uutr_starts) > 0:

        mexon_start[len(uutr_stops) - 1] = uutr_stops[-1]

        mexon_start = mexon_start[len(uutr_stops) - 1:]
        mexon_end = mexon_end[len(uutr_stops) - 1:]

    elif len(dutr_starts) > 0:
        mexon_end[-len(dutr_starts)] = dutr_starts[0]

        mexon_start = mexon_start[:len(exon_start) - len(dutr_starts) + 1]
        mexon_end = mexon_end[:len(exon_start) - len(dutr_starts) + 1]

    uutr_start_str = ';'.join([str(i) for i in uutr_starts])
    uutr_end_str = ';'.join([str(i) for i in uutr_stops])
    dutr_start_str = ';'.join([str(i) for i in dutr_starts])
    dutr_end_str = ';'.join([str(i) for i in dutr_stops])

    exon_start_modified_str = ';'.join([str(i) for i in mexon_start])
    exon_end_modified_str = ';'.join([str(i) for i in mexon_end])

    region_pos_header = uutr_start_str + '|' + uutr_end_str + '|' + dutr_start_str + '|' + dutr_end_str + '|' + exon_start_modified_str + \
                        '|' + exon_end_modified_str + '|' + str(strand) + '|' + str(chr)

    return region_pos_header, mexon_start, mexon_end


def region_sequence(seq, uutr_starts, uutr_stops, dutr_starts, dutr_stops, mexon_start, mexon_end):
    region_seq = 'I' * len(seq)
    region_seq = 'F' * 1000 + region_seq[1000:len(seq) - 1000] + 'f' * 1000

    for s, e in zip(uutr_starts, uutr_stops):
        region_seq = ('U' * (e - s)).join([region_seq[:s], region_seq[e:]])

    for s, e in zip(dutr_starts, dutr_stops):
        region_seq = ('D' * (e - s)).join([region_seq[:s], region_seq[e:]])

    for s, e in zip(mexon_start, mexon_end):
        region_seq = ('E' * (e - s)).join([region_seq[:s], region_seq[e:]])

    return region_seq


def formatting_file(n):

    modified_human_genes_file = open(REGION_FILE_PATH+"\/region_group_" + str(n) + ".fasta", "w")

    for seq_record in SeqIO.parse(SEQ_FILE_PATH+"\group_" + str(n) + ".fasta", "fasta"):
        if GROUP == 'Vertebrate':

            header = (seq_record.id).split('|')
            print(header[2])
            id = header[2]
            sequence = seq_record.seq
            cod_strand = int(header[11])

            GENE_ID, TRANSCRIPT_ID, FIVE_UTR_START, FIVE_UTR_END, THREE_UTR_START, THREE_UTR_END, EXON_STARTS, EXON_ENDS, TSS, EXON_RANK, CHR = genomic_positions_eukaryote(
                cod_strand, header)
        else:
            header = (seq_record.id).split('|')
            print(header[1])
            if SPECIES_NAME == 'lsativa_eg_gene':
                id = header[1] + '_' + header[2] + '_' + header[3]
                del header[1:4]
                header.insert(1, id)
            id = header[1]
            sequence = seq_record.seq
            cod_strand = int(header[9])
            GENE_ID, TRANSCRIPT_ID, FIVE_UTR_START, FIVE_UTR_END, THREE_UTR_START, THREE_UTR_END, EXON_STARTS, EXON_ENDS, TSS, EXON_RANK, CHR = genomic_positions_plant(cod_strand, header)

        if cod_strand == 1:
            uutr_starts, uutr_stops, dutr_starts, dutr_stops = positive_strand_positions(FIVE_UTR_START, FIVE_UTR_END,
                                                                                         THREE_UTR_START, THREE_UTR_END,
                                                                                         TSS)
        else:
            uutr_starts, uutr_stops, dutr_starts, dutr_stops = negative_strand_positions(FIVE_UTR_START, FIVE_UTR_END,
                                                                                         THREE_UTR_START, THREE_UTR_END,
                                                                                         TSS)
        exon_start, exon_end = exon_rank(EXON_RANK, EXON_STARTS, EXON_ENDS, TSS, cod_strand)

        region_pos_header, mexon_start, mexon_end = region_header(exon_start, exon_end, uutr_starts, uutr_stops,
                                                                  dutr_starts, dutr_stops, cod_strand, CHR)
        region_seq = region_sequence(sequence, uutr_starts, uutr_stops, dutr_starts, dutr_stops, mexon_start, mexon_end)

        modified_human_genes_file.write(">" + id + "|" + region_pos_header + "\n")
        modified_human_genes_file.write(region_seq + "\n")
    modified_human_genes_file.close()

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