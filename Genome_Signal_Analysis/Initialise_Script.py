SPECIES_NAME = 'osativa'
GROUP = 'Plant' #'Eukaryote'

Total_Input_Files = 2
COL = [str(i) for i in range(2000)]

ENSEMBLE_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Ensembles_Data\/" + SPECIES_NAME + ".fasta"

DATA_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\/" + SPECIES_NAME
RESULT_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" + SPECIES_NAME


REGION_FILE_PATH = DATA_DIR + "\output_fasta"
SEQ_FILE_PATH = DATA_DIR + "\input_fasta"

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