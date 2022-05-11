SPECIES_NAME = 'h_sapiens' ###'osativa' or 'h_sapiens
GROUP = 'Eukaryote' #'Eukaryote' or 'Plant'

Total_Input_Files = 1
COL = [str(i) for i in range(2000)]
SEQ_TYPE = 'CpG' ### CG_only_pre_mRNA or CG_masked_pre_mRNA or masked_pre_mRNA (this is for repeats masked) or pre_mRNA or mRNA_TSC or mRNA ( This for mRNA script only for now)
NUCLEOSOMAL = True


UCSC_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\UCSC_Data\/" + SPECIES_NAME + ".fasta"
ENSEMBLE_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Ensembles_Data\/" + SPECIES_NAME + ".fasta"
DATA_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\/" + SPECIES_NAME

CG_MASKED_FILE = DATA_DIR + '\CG_masked_flat.fasta' ### Use the SEQ_FILE_PATH to create this file by awk and sed commands
CG_ONLY_FILE = DATA_DIR + '\CpG_only_flat.fasta'

###Di-nucleotide file
GG_ONLY_FILE = DATA_DIR + '\GG_only_flat.fasta'
AG_ONLY_FILE = DATA_DIR + '\AG_only_flat.fasta'
TG_ONLY_FILE = DATA_DIR + '\TG_only_flat.fasta'
CC_ONLY_FILE = DATA_DIR + '\CC_only_flat.fasta'
GC_ONLY_FILE = DATA_DIR + '\GC_only_flat.fasta'
AC_ONLY_FILE = DATA_DIR + '\AC_only_flat.fasta'
TC_ONLY_FILE = DATA_DIR + '\TC_only_flat.fasta'
CA_ONLY_FILE = DATA_DIR + '\CA_only_flat.fasta'
GA_ONLY_FILE = DATA_DIR + '\GA_only_flat.fasta'
AA_ONLY_FILE = DATA_DIR + '\AA_only_flat.fasta'
TA_ONLY_FILE = DATA_DIR + '\TA_only_flat.fasta'
CT_ONLY_FILE = DATA_DIR + '\CT_only_flat.fasta'
GT_ONLY_FILE = DATA_DIR + '\GT_only_flat.fasta'
AT_ONLY_FILE = DATA_DIR + '\AT_only_flat.fasta'
TT_ONLY_FILE = DATA_DIR + '\TT_only_flat.fasta'


if NUCLEOSOMAL:
    RESULT_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" + SPECIES_NAME + r'\NUCLO\/' + SEQ_TYPE
else:
    RESULT_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" + SPECIES_NAME + r'\/' + SEQ_TYPE

REGION_FILE_PATH = DATA_DIR + "\output_fasta"
SEQ_FILE_PATH = DATA_DIR + "\input_fasta"
UCSC_FILE_GROUP_PATH = DATA_DIR + r"\ucsc_fasta"
UCSC_REGION_GROUP_PATH = DATA_DIR + r"\ucsc_region_fasta"

CG_MASKED_SEQ_PATH = DATA_DIR + r"\CG_masked_input_fasta"
CG_MASKED_REGION_PATH = DATA_DIR + r"\CG_masked_output_fasta" ### This is same as the REGION_FILE_PATH

CG_ONLY_SEQ_PATH = DATA_DIR + r"\CpG_only_input_fasta"
CG_ONLY_REGION_PATH = DATA_DIR + r"\CpG_only_output_fasta" ### This is same as the REGION_FILE_PATH
UCSC_SPLIT_RAW_PATH = DATA_DIR + r"\ucsc_split_raw_data"

##DI nucleotide_folders

GG_ONLY_SEQ_PATH = DATA_DIR + r"\GG_only_input_fasta"
AG_ONLY_SEQ_PATH = DATA_DIR + r"\AG_only_input_fasta"
TG_ONLY_SEQ_PATH = DATA_DIR + r"\TG_only_input_fasta"
CC_ONLY_SEQ_PATH = DATA_DIR + r"\CC_only_input_fasta"
GC_ONLY_SEQ_PATH = DATA_DIR + r"\GC_only_input_fasta"
AC_ONLY_SEQ_PATH = DATA_DIR + r"\AC_only_input_fasta"
TC_ONLY_SEQ_PATH = DATA_DIR + r"\TC_only_input_fasta"
CA_ONLY_SEQ_PATH = DATA_DIR + r"\CA_only_input_fasta"
GA_ONLY_SEQ_PATH = DATA_DIR + r"\GA_only_input_fasta"
AA_ONLY_SEQ_PATH = DATA_DIR + r"\AA_only_input_fasta"
TA_ONLY_SEQ_PATH = DATA_DIR + r"\TA_only_input_fasta"
CT_ONLY_SEQ_PATH = DATA_DIR + r"\CT_only_input_fasta"
GT_ONLY_SEQ_PATH = DATA_DIR + r"\GT_only_input_fasta"
AT_ONLY_SEQ_PATH = DATA_DIR + r"\AT_only_input_fasta"
TT_ONLY_SEQ_PATH = DATA_DIR + r"\TT_only_input_fasta"



##### Results File Paths

ABSOLUTE_GC_FILE = RESULT_DIR + r"\Files\absolute_gc_content.csv"
ELEMENT_DENSITY_FILE = RESULT_DIR + r"\Files\elements_density.csv"
AVG_GC_PER_REGION_FILE = RESULT_DIR + r"\Files\Avg_GC_per_element.csv"
GC_PER_REGION_PER_BP_FILE = RESULT_DIR + r"\Files\element_GC_density_perbasepos.csv"

#### Results Charts File Path

ABSOLUTE_GC_CHART = RESULT_DIR + '\Charts\Chart0_absolute_gc_content.png'
ELEMENT_DENSITY_CHART = RESULT_DIR + '\Charts\Chart1_Element_density.png'
AVG_GC_PER_REGION_CHART = RESULT_DIR + '\Charts\Chart2_Avg_GC_Region.png'
GC_PER_REGION_PER_BP_CHART = RESULT_DIR + '\Charts\Chart3_GC_Region_density.png'


### Derived Charts Location

ABS_GC_DENSITY_BY_REGION_CHART = RESULT_DIR + '\Charts\Chart4_GC_Density_and_element_density.png'
DERIVED_ABSOLUTE_GC_CHART = RESULT_DIR + '\Charts\Chart6_Derived_Absolute_GC.png'
AVG_GC_DENSITY_BY_REGION_CHART = RESULT_DIR + '\Charts\Chart5_Avg_GC_Density_and_element_density.png'
DERIVED_AVG_GC_CHART = RESULT_DIR + '\Charts\Chart7_Derived_Avg_Absolute_GC.png'

CENTERED_GC_PER_REGION_PER_BP_CHART= RESULT_DIR + '\Charts\Chart8_Centered_GC_Region_density.png'
CENTERED_GC_DENSITY_BY_REGION_CHART= RESULT_DIR + '\Charts\Chart9_Centered_GC_density_and element_density.png'
DERIVED_CENTERED_GC_CHART = RESULT_DIR + '\Charts\Chart10_Centered_Derived_Absolute_GC.png'

COMPARISON_DERIVED_ABS_GC = RESULT_DIR + '\Charts\Chart11_Derived_and_Absolute_GC.png'

### Temporary Files Location

AVG_GC_DENSITY_BY_REGION_FILE = RESULT_DIR + '\Files\Temp_Chart_5.csv'
DERIVED_AVG_GC_FILE = RESULT_DIR + '\Files\Temp_Chart_7.csv'
CENTERED_GC_PER_REGION_PER_BP_FILE = RESULT_DIR+'\Files\Temp_Chart_8.csv'
CENTERED_GC_DENSITY_BY_REGION_FILE = RESULT_DIR+'\Files\Temp_Chart_9.csv'
DERIVED_CENTERED_GC_FILE = RESULT_DIR+'\Files\Temp_Chart_10.csv'

#### Combine Charts Location

ABSOLUTE_GC_AND_REGION_DENSITY_CHARTS = RESULT_DIR+'\Charts\combine.png'
INTRA_AND_INTER_REGION_CHARTS = RESULT_DIR+'\Charts\combine_1.png'

###CPG CHART LOCATION

ABSOLUTE_CPG_CHART = RESULT_DIR + '\Charts\Chart0_absolute_CPG_content.png'
REGION_DENSITY_CHART = RESULT_DIR + '\Charts\Chart1_region_density.png'
AVG_CPG_PER_REGION_CHART = RESULT_DIR + '\Charts\Chart2_Avg_CPG_Region.png'
CPG_PER_REGION_PER_BP_CHART = RESULT_DIR + '\Charts\Chart3_CPG_Region_density.png'


###CPG FILE LOCATION

ABSOLUTE_CPG_FILE = RESULT_DIR + '\Charts\Chart0_absolute_CPG_content.csv'
REGION_DENSITY_FILE = RESULT_DIR + '\Charts\Chart1_region_density.csv'
AVG_CPG_PER_REGION_FILE = RESULT_DIR + '\Charts\Chart2_Avg_CPG_Region.csv'
CPG_PER_REGION_PER_BP_FILE = RESULT_DIR + '\Charts\Chart3_CPG_Region_density.csv'



####MRNA FILE NAMES####################################################################
##### Results File Paths

# ABSOLUTE_GC_FILE = RESULT_DIR + r"\Files\MRNA_absolute_gc_content.csv"
# ELEMENT_DENSITY_FILE = RESULT_DIR + r"\Files\MRNA_elements_density.csv"
# AVG_GC_PER_REGION_FILE = RESULT_DIR + r"\Files\MRNA_Avg_GC_per_element.csv"
# GC_PER_REGION_PER_BP_FILE = RESULT_DIR + r"\Files\MRNA_element_GC_density_perbasepos.csv"
#
# #### Results Charts File Path
#
# ABSOLUTE_GC_CHART = RESULT_DIR + '\Charts\MRNA_Chart0_absolute_gc_content.png'
# ELEMENT_DENSITY_CHART = RESULT_DIR + '\Charts\MRNA_Chart1_Element_density.png'
# AVG_GC_PER_REGION_CHART = RESULT_DIR + '\Charts\MRNA_Chart2_Avg_GC_Region.png'
# GC_PER_REGION_PER_BP_CHART = RESULT_DIR + '\Charts\MNRA_Chart3_GC_Region_density.png'
#
# ### Derived Charts Location
#
# ABS_GC_DENSITY_BY_REGION_CHART = RESULT_DIR + '\Charts\MRNA_Chart4_GC_Density_and_element_density.png'
# DERIVED_ABSOLUTE_GC_CHART = RESULT_DIR + '\Charts\MRNA_Chart6_Derived_Absolute_GC.png'
# AVG_GC_DENSITY_BY_REGION_CHART = RESULT_DIR + '\Charts\MRNA_Chart5_Avg_GC_Density_and_element_density.png'
# DERIVED_AVG_GC_CHART = RESULT_DIR + '\Charts\MRNA_Chart7_Derived_Avg_Absolute_GC.png'
#
# CENTERED_GC_PER_REGION_PER_BP_CHART= RESULT_DIR + '\Charts\MRNA_Chart8_Centered_GC_Region_density.png'
# CENTERED_GC_DENSITY_BY_REGION_CHART= RESULT_DIR + '\Charts\MRNA_Chart9_Centered_GC_density_and element_density.png'
# DERIVED_CENTERED_GC_CHART = RESULT_DIR + '\Charts\MRNA_Chart10_Centered_Derived_Absolute_GC.png'
#
# COMPARISON_DERIVED_ABS_GC = RESULT_DIR + '\Charts\MRNA_Chart11_Derived_and_Absolute_GC.png'
#
# ### Temporary Files Location
#
# AVG_GC_DENSITY_BY_REGION_FILE = RESULT_DIR + '\Files\MRNA_Temp_Chart_5.csv'
# DERIVED_AVG_GC_FILE = RESULT_DIR + '\Files\MRNA_Temp_Chart_7.csv'
# CENTERED_GC_PER_REGION_PER_BP_FILE = RESULT_DIR+'\Files\MRNA_Temp_Chart_8.csv'
# CENTERED_GC_DENSITY_BY_REGION_FILE = RESULT_DIR+'\Files\MRNA_Temp_Chart_9.csv'
# DERIVED_CENTERED_GC_FILE = RESULT_DIR+'\Files\MRNA_Temp_Chart_10.csv'
#
# #### Combine Charts Location
#
# ABSOLUTE_GC_AND_REGION_DENSITY_CHARTS = RESULT_DIR+'\Charts\MRNA_combine.png'
# INTRA_AND_INTER_REGION_CHARTS = RESULT_DIR+'\Charts\MRNA_combine_1.png'