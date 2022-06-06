import sys
print(sys.argv)

SPECIES_NAME = sys.argv[1]
# SPECIES_NAME = 'ggallus_gene_ensembl' ###'osativa' or 'h_sapiens (Can be read from Arguments)
GROUP = 'Vertebrate' #'Vertebrate' or 'Plant' or 'Yeast'

REGION_DICT = {'U': 1, 'D': 2, 'E': 3, 'I': 4, 'F': 5, 'f': 6}
OPP_REGION_DICT = {v: k for k, v in REGION_DICT.items()}

COL = [str(i) for i in range(2000)]
SEQ_TYPE = 'mRNA' # pre_mRNA, mRNA, mRNA_ATG

ENSEMBLE_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Ensembles_Data\/" + SPECIES_NAME + ".fasta"
DATA_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\/" + SPECIES_NAME

RESULT_NUCLO_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/"+GROUP+ '\/' + SPECIES_NAME + r'\NUCLO\/' + SEQ_TYPE
RESULT_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" + GROUP + r"\/" + SPECIES_NAME + r"\/" + SEQ_TYPE

if SEQ_TYPE == 'pre_mRNA':

    SEQ_FILE_PATH = DATA_DIR + r"\seq_dir"
    REGION_FILE_PATH = DATA_DIR + r"\region_dir"
    Total_Input_Files = 2
    # Total_Input_Files = len([name for name in os.listdir(SEQ_FILE_PATH) if os.path.isfile(os.path.join(SEQ_FILE_PATH, name))])

if SEQ_TYPE in ['mRNA', 'mRNA_ATG']:

    pre_mrna_SEQ_FILE = DATA_DIR + r"\seq_dir"
    pre_mrna_REGION_FILE = DATA_DIR + r"\region_dir"

    SEQ_FILE_PATH = DATA_DIR + r"\mrna_seq_dir"
    REGION_FILE_PATH = DATA_DIR + r"\mrna_region_dir"

    Total_Input_Files = 2
    # Total_Input_Files = len([name for name in os.listdir(pre_mrna_SEQ_FILE) if os.path.isfile(os.path.join(pre_mrna_SEQ_FILE, name))])


##### Results File Paths

ABSOLUTE_GC_FILE = RESULT_DIR + r"\Files\absolute_gc_content.csv" ##keep
NUCLO_ABSOLUTE_GC_FILE = RESULT_NUCLO_DIR + r"\Files\absolute_gc_content.csv" ##keep


ELEMENT_DENSITY_FILE = RESULT_DIR + r"\Files\elements_density.csv" ##keep
NUCLO_ELEMENT_DENSITY_FILE = RESULT_NUCLO_DIR + r"\Files\elements_density.csv" ##keep


AVG_GC_PER_REGION_FILE = RESULT_DIR + r"\Files\Avg_GC_per_element.csv" ##keep
INTER_REGION_SIGNAL_DIST_FILE = RESULT_DIR + '\Files\inter_region_signal_per_region.csv' ###keep
INTER_REGION_SIGNAL_FILE = RESULT_DIR + '\Files\inter_region_signal.csv' ##keep
NUCLO_AVG_GC_PER_REGION_FILE = RESULT_NUCLO_DIR + r"\Files\Avg_GC_per_element.csv" ##keep
NUCLO_INTER_REGION_SIGNAL_DIST_FILE = RESULT_NUCLO_DIR + '\Files\inter_region_signal_per_region.csv' ##keep
NUCLO_INTER_REGION_SIGNAL_FILE = RESULT_NUCLO_DIR + '\Files\inter_region_signal.csv' ##keep


NORM_INTRA_REGION_SIGNAL_DIST_FILE = RESULT_DIR + r"\Files\norm_intra_region_signal_per_region.csv" ##keep
CENTERED_NORM_INTRA_REGION_SIGNAL_DIST_FILE = RESULT_DIR+'\Files\Centered_norm_intra_region_signal_per_region.csv' ##keep
INTRA_REGION_SIGNAL_DIST_FILE = RESULT_DIR + '\Files\intra_region_signal_per_region.csv' ##keep
CENTERED_INTRA_REGION_SIGNAL_DIST_FILE = RESULT_DIR+'\Files\Centered_intra_region_signal_per_region.csv' ##keep
INTRA_REGION_SIGNAL_FILE = RESULT_DIR+'\Files\intra_region_signal.csv' ##keep

NUCLO_NORM_INTRA_REGION_SIGNAL_DIST_FILE = RESULT_NUCLO_DIR + r"\Files\norm_intra_region_signal_per_region.csv" ##keep
NUCLO_CENTERED_NORM_INTRA_REGION_SIGNAL_DIST_FILE = RESULT_NUCLO_DIR+'\Files\Centered_norm_intra_region_signal_per_region.csv' ##keep
NUCLO_INTRA_REGION_SIGNAL_DIST_FILE = RESULT_NUCLO_DIR + '\Files\intra_region_signal_per_region.csv' ##keep
NUCLO_CENTERED_INTRA_REGION_SIGNAL_DIST_FILE = RESULT_NUCLO_DIR+'\Files\Centered_intra_region_signal_per_region.csv' ##keep
NUCLO_INTRA_REGION_SIGNAL_FILE = RESULT_NUCLO_DIR+'\Files\intra_region_signal.csv' ##keep

#### Results Charts File Path

ABSOLUTE_GC_CHART = RESULT_DIR + '\SEL_CHARTS\ABSOLUTE_GC_CONTENT.png' ##keep
NUCLO_ABSOLUTE_GC_CHART = RESULT_NUCLO_DIR + '\SEL_CHARTS\ABSOLUTE_GC_CONTENT.png' ##keep

ABSOLUTE_and_SIGNAL_GC_LINE_CHART = RESULT_DIR+'\SEL_CHARTS\ABSOLUTE_GC_AND_DERIVED_SIGNAL.png' ##keep
ELEMENT_DENSITY_LINE_CHART = RESULT_DIR+'\SEL_CHARTS\ELEMENT_DENSITY.png' ##keep

AVG_GC_REGION_CHART = RESULT_DIR+'\SEL_CHARTS\AVG_GC_REGION.png' ##keep
INTER_SIGNAL_REGION_DIST_CHART = RESULT_DIR+'\SEL_CHARTS\INTER_REGION_SIGNAL_DIST.png' ##keep

NORM_INTRA_SIGNAL_REGION_DIST_CHART = RESULT_DIR+'\SEL_CHARTS\/NORM_INTRA_REGION_SIGNAL_DIST.png' ##keep
INTRA_SIGNAL_REGION_DIST_CHART = RESULT_DIR+'\SEL_CHARTS\INTRA_REGION_SIGNAL_DIST.png' ##keep

NUCLO_ABSOLUTE_and_SIGNAL_GC_LINE_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\ABSOLUTE_GC_AND_DERIVED_SIGNAL.png' ##keep
NUCLO_ELEMENT_DENSITY_LINE_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\ELEMENT_DENSITY.png' ##keep

NUCLO_AVG_GC_REGION_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\AVG_GC_REGION.png' ##keep
NUCLO_INTER_SIGNAL_REGION_DIST_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\INTER_REGION_SIGNAL_DIST.png'##keep

NUCLO_NORM_INTRA_SIGNAL_REGION_DIST_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\/NORM_INTRA_REGION_SIGNAL_DIST.png' ##keep
NUCLO_INTRA_SIGNAL_REGION_DIST_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\INTRA_REGION_SIGNAL_DIST.png' ##keep
