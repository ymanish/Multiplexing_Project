import sys
import os
print(sys.argv)

AA_DICT= {'A': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7,
           'I': 8, 'K': 9, 'L': 10, 'M': 11, 'N': 12, 'P': 13, 'Q': 14,
           'R': 15, 'S': 16, 'T': 17, 'U': 18, 'V': 19, 'W': 20,
           'Y': 21, 'X': 22, 'B': 23, 'Z': 24, 'J': 25, 'O': 26,
           'x': 27, 'u': 28, 'i': 29, 'd': 30, 'f': 31, 'z': 32}
OPP_AA_DICT = {v: k for k, v in AA_DICT.items()}

SPECIES_NAME = sys.argv[1]
print(SPECIES_NAME)
# SPECIES_NAME = 'h_sapiens' ###'osativa' or 'h_sapiens (Can be read from Arguments)
GROUP = 'Vertebrate' #'Vertebrate' or 'Plant' or 'Yeast'

COL = [str(i) for i in range(2000)]
SEQ_TYPE = sys.argv[2] ## Amino_Acid, Amino_Acid_mRNA, Amino_Acid_mRNA_ATG

ENSEMBLE_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Ensembles_Data\/" + SPECIES_NAME + ".fasta"
DATA_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\/" + SPECIES_NAME

RESULT_NUCLO_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/"+GROUP+ '\/' + SPECIES_NAME + r'\NUCLO\/' + SEQ_TYPE
RESULT_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" +GROUP+ '\/' +SPECIES_NAME + r'\/' + SEQ_TYPE


pre_mrna_SEQ_FILE_PATH = DATA_DIR + r"\seq_dir"
pre_mrna_REGION_FILE_PATH = DATA_DIR + r"\region_dir"
AA_FILE_PATH = DATA_DIR + r"\aa_dir"

# Total_Input_Files = len([name for name in os.listdir(pre_mrna_SEQ_FILE_PATH) if os.path.isfile(os.path.join(pre_mrna_SEQ_FILE_PATH, name))])
Total_Input_Files = 2

mrna_SEQ_FILE_PATH = DATA_DIR + r"\mrna_seq_dir"
mrna_REGION_FILE_PATH = DATA_DIR + r"\mrna_region_dir"
mrna_AA_FILE_PATH = DATA_DIR + r"\mrna_aa_dir"

if SEQ_TYPE == 'Amino_Acid':

    CDS_INTER_REGION_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" +GROUP+ '\/' +SPECIES_NAME + r'\pre_mRNA\Files\inter_region_signal_per_region.csv'
    NUCLO_CDS_INTER_REGION_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" +GROUP+ '\/' +SPECIES_NAME + r'\NUCLO\pre_mRNA\Files\inter_region_signal_per_region.csv'

if SEQ_TYPE == 'Amino_Acid_mRNA':
    CDS_INTER_REGION_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" +GROUP+ '\/' +SPECIES_NAME + r'\mRNA\Files\inter_region_signal_per_region.csv'
    NUCLO_CDS_INTER_REGION_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" +GROUP+ '\/' +SPECIES_NAME + r'\NUCLO\mRNA\Files\inter_region_signal_per_region.csv'

if SEQ_TYPE == 'Amino_Acid_mRNA_ATG':
    CDS_INTER_REGION_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" + GROUP + '\/' + SPECIES_NAME + r'\mRNA_ATG\Files\inter_region_signal_per_region.csv'
    NUCLO_CDS_INTER_REGION_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" + GROUP + '\/' + SPECIES_NAME + r'\NUCLO\mRNA_ATG\Files\inter_region_signal_per_region.csv'


##### Results File Paths

'AMINO ACID CHOICE SIGNAL FILES'

AMINO_ACID_DENSITY_FILE = RESULT_DIR + r"\Files\aa_density.csv"
AVG_GC_AA_FILE = RESULT_DIR + r"\Files\avg_gc_aa.csv"
AA_CHOICE_SIGNAL_FILE = RESULT_DIR + r"\Files\aa_choice_signal_dist.csv"
TOTAL_AA_CHOICE_SIGNAL_FILE = RESULT_DIR + r"\Files\total_aa_choice_signal.csv"
CORRECTED_TOTAL_AA_CHOICE_SIGNAL_FILE = RESULT_DIR + r"\Files\corrected_total_aa_choice_signal.csv"

NUCLO_AMINO_ACID_DENSITY_FILE = RESULT_NUCLO_DIR + r"\Files\aa_density.csv"
NUCLO_AVG_GC_AA_FILE = RESULT_NUCLO_DIR + r"\Files\avg_gc_aa.csv"
NUCLO_AA_CHOICE_SIGNAL_FILE = RESULT_NUCLO_DIR + r"\Files\aa_choice_signal_dist.csv"
NUCLO_TOTAL_AA_CHOICE_SIGNAL_FILE = RESULT_NUCLO_DIR + r"\Files\total_aa_choice_signal.csv"
NUCLO_CORRECTED_TOTAL_AA_CHOICE_SIGNAL_FILE = RESULT_NUCLO_DIR + r"\Files\corrected_total_aa_choice_signal.csv"


'SYNONYMOUS CODON SIGNAL FILES'
SYN_CODON_SIGNAL_FILE = RESULT_DIR + r"\Files\syn_codon_signal_dist.csv"
CEN_SYN_CODON_SIGNAL_FILE =  RESULT_DIR + r"\Files\centered_syn_codon_signal_dist.csv"
NORMALISED_SYN_CODON_SIGNAL_FILE= RESULT_DIR + r"\Files\normalised_syn_codon_signal_dist.csv"
CENT_NORMALISED_SYN_CODON_SIGNAL_FILE= RESULT_DIR + r"\Files\centered_normalised_syn_codon_signal_dist.csv"
TOTAL_SYN_CODON_SIGNAL_FILE = RESULT_DIR + r"\Files\total_syn_signal.csv"

NUCLO_SYN_CODON_SIGNAL_FILE = RESULT_NUCLO_DIR + r"\Files\syn_codon_signal_dist.csv"
NUCLO_CEN_SYN_CODON_SIGNAL_FILE =  RESULT_NUCLO_DIR + r"\Files\centered_syn_codon_signal_dist.csv"
NUCLO_NORMALISED_SYN_CODON_SIGNAL_FILE= RESULT_NUCLO_DIR + r"\Files\normalised_syn_codon_signal_dist.csv"
NUCLO_CENT_NORMALISED_SYN_CODON_SIGNAL_FILE= RESULT_NUCLO_DIR + r"\Files\centered_normalised_syn_codon_signal_dist.csv"
NUCLO_TOTAL_SYN_CODON_SIGNAL_FILE = RESULT_NUCLO_DIR + r"\Files\total_syn_signal.csv"

#### Results Charts File Path
'AMINO ACID CHOICE SIGNAL CHARTS'

AA_CHOICE_SIGNAL_CHART = RESULT_DIR + r"\Charts\aa_choice_signal_dist.png"
TOTAL_AA_CHOICE_SIGNAL_CHART = RESULT_DIR + r"\Charts\total_aa_choice_signal.png"
CORRECTED_TOTAL_AA_CHOICE_SIGNAL_CHART=  RESULT_DIR + r"\Charts\corrected_total_aa_choice_signal.png"


NUCLO_AA_CHOICE_SIGNAL_CHART = RESULT_NUCLO_DIR + r"\Charts\aa_choice_signal_dist.png"
NUCLO_TOTAL_AA_CHOICE_SIGNAL_CHART = RESULT_NUCLO_DIR + r"\Charts\total_aa_choice_signal.png"
NUCLO_CORRECTED_TOTAL_AA_CHOICE_SIGNAL_CHART=  RESULT_NUCLO_DIR + r"\Charts\corrected_total_aa_choice_signal.png"


'SYNONYMOUS CODON SIGNAL CHARTS'
CEN_SYN_CODON_SIGNAL_CHART = RESULT_DIR + r"\Charts\centered_syn_codon_signal_dist.png"
CENT_NORMALISED_SYN_CODON_SIGNAL_CHART= RESULT_DIR + r"\Charts\centered_normalised_syn_codon_signal_dist.png"
TOTAL_SYN_CODON_SIGNAL_CHART = RESULT_DIR + r"\Charts\total_syn_signal.png"

NUCLO_CEN_SYN_CODON_SIGNAL_CHART = RESULT_NUCLO_DIR + r"\Charts\centered_syn_codon_signal_dist.png"
NUCLO_CENT_NORMALISED_SYN_CODON_SIGNAL_CHART= RESULT_NUCLO_DIR + r"\Charts\centered_normalised_syn_codon_signal_dist.png"
NUCLO_TOTAL_SYN_CODON_SIGNAL_CHART = RESULT_NUCLO_DIR + r"\Charts\total_syn_signal.png"

'SELECTED CHARTS OF AA SIGNAL PROFILING'
AMINO_ACID_DENSITY_CHART = RESULT_DIR + r"\SEL_CHARTS\aa_density.png"
AVG_GC_AA_CHART =  RESULT_DIR + r"\SEL_CHARTS\avg_gc_aa.png"
AA_CHOICE_AND_SYN_SIGNAL = RESULT_DIR + r"\SEL_CHARTS\aa_choice_and_syn_signal.png"



NUCLO_AMINO_ACID_DENSITY_CHART = RESULT_NUCLO_DIR + r"\SEL_CHARTS\aa_density.png"
NUCLO_AVG_GC_AA_CHART = RESULT_NUCLO_DIR + r"\SEL_CHARTS\avg_gc_aa.png"
NUCLO_AA_CHOICE_AND_SYN_SIGNAL = RESULT_NUCLO_DIR + r"\SEL_CHARTS\aa_choice_and_syn_signal.png"


