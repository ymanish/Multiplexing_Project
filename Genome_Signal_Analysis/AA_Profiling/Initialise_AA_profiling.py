import os.path
import sys

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
SEQ_TYPE = 'Amino_Acid'

ENSEMBLE_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Ensembles_Data\/" + SPECIES_NAME + ".fasta"
DATA_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\/" + SPECIES_NAME

RESULT_NUCLO_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/"+GROUP+ '\/' + SPECIES_NAME + r'\NUCLO\/' + SEQ_TYPE
RESULT_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" +GROUP+ '\/' +SPECIES_NAME + r'\/' + SEQ_TYPE

SEQ_FILE_PATH = DATA_DIR + r"\seq_dir"
REGION_FILE_PATH = DATA_DIR + r"\region_dir"
AA_FILE_PATH = DATA_DIR + r"\aa_dir"

# Total_Input_Files = len([name for name in os.listdir(SEQ_FILE_PATH) if os.path.isfile(os.path.join(SEQ_FILE_PATH, name))])
Total_Input_Files = 15
##### Results File Paths

AMINO_ACID_DENSITY_FILE = RESULT_DIR + r"\Files\aa_density.csv"
AVG_GC_AA_FILE = RESULT_DIR + r"\Files\avg_gc_aa.csv"
AA_CHOICE_SIGNAL_FILE = RESULT_DIR + r"\Files\aa_choice_signal.csv"
TOTAL_INTER_SIGNAL_FILE = RESULT_DIR + r"\Files\inter_signal.csv"
SYN_CODON_SIGNAL_FILE = RESULT_DIR + r"\Files\syn_codon_signal.csv"
CEN_SYN_CODON_SIGNAL_FILE =  RESULT_DIR + r"\Files\centered_syn_codon_signal.csv"
NORMALISED_SYN_CODON_SIGNAL_FILE= RESULT_DIR + r"\Files\normalised_syn_codon_signal.csv"
CENT_NORMALISED_SYN_CODON_SIGNAL_FILE= RESULT_DIR + r"\Files\centered_normalised_syn_codon_signal.csv"
TOTAL_INTRA_SIGNAL_FILE = RESULT_DIR + r"\Files\intra_signal.csv"

#### Results Charts File Path

AMINO_ACID_DENSITY_CHART = RESULT_DIR + r"\Charts\aa_density.png"
AVG_GC_AA_CHART =  RESULT_DIR + r"\Charts\avg_gc_aa.png"
AA_CHOICE_SIGNAL_CHART = RESULT_DIR + r"\Charts\aa_choice_signal.png"
TOTAL_INTER_SIGNAL_CHART = RESULT_DIR + r"\Charts\inter_signal.png"
SYN_CODON_SIGNAL_CHART = RESULT_DIR + r"\Charts\syn_codon_signal.png"
CEN_SYN_CODON_SIGNAL_CHART = RESULT_DIR + r"\Charts\centered_syn_codon_signal.png"
NORMALISED_SYN_CODON_SIGNAL_CHART = RESULT_DIR + r"\Charts\normalised_syn_codon_signal.png"
CENT_NORMALISED_SYN_CODON_SIGNAL_CHART= RESULT_DIR + r"\Charts\centered_normalised_syn_codon_signal.png"
TOTAL_INTRA_SIGNAL_CHART = RESULT_DIR + r"\Charts\intra_signal.png"
