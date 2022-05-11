import os.path
import sys

print(sys.argv)

SPECIES_NAME = sys.argv[1]
# SPECIES_NAME = 'h_sapiens' ###'osativa' or 'h_sapiens (Can be read from Arguments)
GROUP = 'Vertebrate' #'Vertebrate' or 'Plant' or 'Yeast'

COL = [str(i) for i in range(2000)]
SEQ_TYPE = 'pre_mRNA'

ENSEMBLE_FILE = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Ensembles_Data\/" + SPECIES_NAME + ".fasta"
DATA_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Data\/" + SPECIES_NAME

RESULT_NUCLO_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/"+GROUP+ '\/' + SPECIES_NAME + r'\NUCLO\/' + SEQ_TYPE
RESULT_DIR = r"C:\Users\maya620d\PycharmProjects\Multiplexing\Results\/" + GROUP + r"\/" + SPECIES_NAME + r"\/" + SEQ_TYPE

SEQ_FILE_PATH = DATA_DIR + r"\seq_dir"
REGION_FILE_PATH = DATA_DIR + r"\region_dir"
Total_Input_Files = 15
# Total_Input_Files = len([name for name in os.listdir(SEQ_FILE_PATH) if os.path.isfile(os.path.join(SEQ_FILE_PATH, name))])

##### Results File Paths

ABSOLUTE_GC_FILE = RESULT_DIR + r"\Files\absolute_gc_content.csv"
ELEMENT_DENSITY_FILE = RESULT_DIR + r"\Files\elements_density.csv"
AVG_GC_PER_REGION_FILE = RESULT_DIR + r"\Files\Avg_GC_per_element.csv"
GC_PER_REGION_PER_BP_FILE = RESULT_DIR + r"\Files\element_GC_density_perbasepos.csv"
GC_PER_REGION_PER_BP_ALL_TRANS_FILE = RESULT_DIR + '\Files\GC_Density_and_element_density.csv'
AVG_GC_DENSITY_BY_REGION_FILE = RESULT_DIR + '\Files\Avg_GC_Density_and_element_density.csv'
DERIVED_CUM_AVG_GC_FILE = RESULT_DIR + '\Files\Derived_Avg_Absolute_GC.csv'
DERIVED_ABSOLUTE_GC_FILE = RESULT_DIR + '\Files\Derived_Absolute_GC.csv'


NUCLO_ABSOLUTE_GC_FILE = RESULT_NUCLO_DIR + r"\Files\absolute_gc_content.csv"
NUCLO_ELEMENT_DENSITY_FILE = RESULT_NUCLO_DIR + r"\Files\elements_density.csv"
NUCLO_GC_PER_REGION_PER_BP_FILE = RESULT_NUCLO_DIR + r"\Files\element_GC_density_perbasepos.csv"
NUCLO_GC_PER_REGION_PER_BP_ALL_FILE = RESULT_NUCLO_DIR + '\Files\GC_Density_and_element_density.csv'
NUCLO_AVG_GC_DENSITY_BY_REGION_FILE = RESULT_NUCLO_DIR + '\Files\Avg_GC_Density_and_element_density.csv'
NUCLO_DERIVED_CUM_AVG_GC_FILE = RESULT_NUCLO_DIR + '\Files\Derived_Avg_Absolute_GC.csv'
NUCLO_DERIVED_ABSOLUTE_GC_FILE = RESULT_NUCLO_DIR + '\Files\Derived_Absolute_GC.csv'


#### Results Charts File Path

ABSOLUTE_GC_CHART = RESULT_DIR + '\Charts\Chart0_absolute_gc_content.png'
ELEMENT_DENSITY_CHART = RESULT_DIR + '\Charts\Chart1_Element_density.png'
AVG_GC_PER_REGION_CHART = RESULT_DIR + '\Charts\Chart2_Avg_GC_Region.png'
GC_PER_REGION_PER_BP_CHART = RESULT_DIR + '\Charts\Chart3_GC_Region_density.png'
GC_PER_REGION_PER_BP_ALL_TRANS_CHART = RESULT_DIR + '\Charts\Chart4_GC_Density_and_element_density.png'
AVG_GC_DENSITY_BY_REGION_CHART = RESULT_DIR + '\Charts\Chart5_Avg_GC_Density_and_element_density.png'
DERIVED_CUM_AVG_GC_CHART = RESULT_DIR + '\Charts\Chart7_Derived_Avg_Absolute_GC.png'
DERIVED_ABSOLUTE_GC_CHART = RESULT_DIR + '\Charts\Chart6_Derived_Absolute_GC.png'


NUCLO_ABSOLUTE_GC_CHART = RESULT_NUCLO_DIR + '\Charts\Chart0_absolute_gc_content.png'
NUCLO_ELEMENT_DENSITY_CHART = RESULT_NUCLO_DIR + '\Charts\Chart1_Element_density.png'
NUCLO_GC_PER_REGION_PER_BP_CHART = RESULT_NUCLO_DIR + '\Charts\Chart3_GC_Region_density.png'
NUCLO_GC_PER_REGION_PER_BP_ALL_CHART = RESULT_NUCLO_DIR + '\Charts\Chart4_GC_Density_and_element_density.png'
NUCLO_AVG_GC_DENSITY_BY_REGION_CHART = RESULT_NUCLO_DIR + '\Charts\Chart5_Avg_GC_Density_and_element_density.png'
NUCLO_DERIVED_CUM_AVG_GC_CHART = RESULT_NUCLO_DIR + '\Charts\Chart7_Derived_Avg_Absolute_GC.png'
NUCLO_DERIVED_ABSOLUTE_GC_CHART = RESULT_NUCLO_DIR + '\Charts\Chart6_Derived_Absolute_GC.png'


### Derived Charts Location

CENTERED_GC_PER_REGION_PER_BP_CHART = RESULT_DIR + '\Charts\Chart8_Centered_GC_Region_density.png'
CENTERED_GC_DENSITY_BY_REGION_CHART = RESULT_DIR + '\Charts\Chart9_Centered_GC_density_and element_density.png'
DERIVED_CENTERED_GC_CHART = RESULT_DIR + '\Charts\Chart10_Centered_Derived_Absolute_GC.png'

NUCLO_CENTERED_GC_PER_REGION_PER_BP_CHART = RESULT_NUCLO_DIR + '\Charts\Chart8_Centered_GC_Region_density.png'
NUCLO_CENTERED_GC_DENSITY_BY_REGION_CHART = RESULT_NUCLO_DIR + '\Charts\Chart9_Centered_GC_density_and element_density.png'
NUCLO_DERIVED_CENTERED_GC_CHART = RESULT_NUCLO_DIR + '\Charts\Chart10_Centered_Derived_Absolute_GC.png'

### Derived Files Location

CENTERED_GC_PER_REGION_PER_BP_FILE = RESULT_DIR+'\Files\Centered_GC_Region_density.csv'
CENTERED_GC_DENSITY_BY_REGION_FILE = RESULT_DIR+'\Files\Centered_GC_density_and element_density.csv'
DERIVED_CENTERED_GC_FILE = RESULT_DIR+'\Files\Centered_Derived_Absolute_GC.csv'

NUCLO_CENTERED_GC_PER_REGION_PER_BP_FILE = RESULT_NUCLO_DIR+'\Files\Centered_GC_Region_density.csv'
NUCLO_CENTERED_GC_DENSITY_BY_REGION_FILE = RESULT_NUCLO_DIR+'\Files\Centered_GC_density_and element_density.csv'
NUCLO_DERIVED_CENTERED_GC_FILE = RESULT_NUCLO_DIR+'\Files\Centered_Derived_Absolute_GC.csv'

#### FINAL GC PROFILING CHARTS LOCATION

ABSOLUTE_GC_LINE_CHART = RESULT_DIR+'\SEL_CHARTS\ABSOLUTE_GC_CONTENT.png'
ELEMENT_DENSITY_LINE_CHART = RESULT_DIR+'\SEL_CHARTS\ELEMENT_DENSITY.png'
DERIVED_CUM_GC_LINE_CHART = RESULT_DIR+'\SEL_CHARTS\DERIVED_CUM_GC_CONTENT.png'
# INTRA_REGION_CHARTS = RESULT_DIR+'\SEL_CHARTS\INTRA_REGIONAL_CHARTS.png'
# INTER_REGION_CHARTS = RESULT_DIR+'\SEL_CHARTS\INTER_REGIONAL_CHARTS.png'

INTRA_SIGNAL_REGION_DIST_PERBP_CHART = RESULT_DIR+'\SEL_CHARTS\INTRA_SIGNAL_REGION_DIST_PERBP.png'
INTRA_SIGNAL_REGION_DIST_CHART = RESULT_DIR+'\SEL_CHARTS\INTRA_SIGNAL_REGION_DIST.png'
AVG_GC_REGION_CHART = RESULT_DIR+'\SEL_CHARTS\AVG_GC_REGION.png'
INTER_SIGNAL_REGION_DIST_CHART = RESULT_DIR+'\SEL_CHARTS\INTER_SIGNAL_REGION_DIST.png'


NUCLO_ABSOLUTE_GC_LINE_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\ABSOLUTE_GC_CONTENT.png'
NUCLO_ELEMENT_DENSITY_LINE_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\ELEMENT_DENSITY.png'
NUCLO_DERIVED_CUM_GC_LINE_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\DERIVED_CUM_GC_CONTENT.png'
# NUCLO_INTRA_REGION_CHARTS = RESULT_NUCLO_DIR+'\SEL_CHARTS\INTRA_REGIONAL_CHARTS.png'
# NUCLO_INTER_REGION_CHARTS = RESULT_NUCLO_DIR+'\SEL_CHARTS\INTER_REGIONAL_CHARTS.png'


NUCLO_INTRA_SIGNAL_REGION_DIST_PERBP_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\INTRA_SIGNAL_REGION_DIST_PERBP.png'
NUCLO_INTRA_SIGNAL_REGION_DIST_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\INTRA_SIGNAL_REGION_DIST.png'
NUCLO_AVG_GC_REGION_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\AVG_GC_REGION.png'
NUCLO_INTER_SIGNAL_REGION_DIST_CHART = RESULT_NUCLO_DIR+'\SEL_CHARTS\INTER_SIGNAL_REGION_DIST.png'