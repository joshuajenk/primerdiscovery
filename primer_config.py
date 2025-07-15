# Thermodynamic thresholds
Def_Melt_Tm_LOWER = 45 #degrees Celsius
Def_Melt_Tm_UPPER = 65 #degrees Celsius
Def_GC_LOWER = 40 #GC%
Def_GC_UPPER = 60 #GC%
Def_Hairpin_Min = 3 #integer
Def_SelfDimer_Min = 5 #integer

# Hydrolysis probe design constraints
Def_HP_Length = 20 #bp
Def_HP_Tm_LOWER = 58 #degrees Celsius
Def_HP_Tm_UPPER = 70 #degrees Celsius
Def_HP_GC_LOWER = 40 #GC%
Def_HP_GC_UPPER = 60 #GC%
#HYDROLYSIS PROBE - FLUORO = 'FAM'
#HYDROLYSIS PROBE - QUENCHER = 'BHQ1'

# Primer Thresholds
Primer3_PRODUCT_SIZE_RANGE = [(100,500)] #bp: keep umder 1000 bp
Primer3_PRIMER_SIZE_RANGE   = (12,30) #bp: keep between 10 - 30 bp
Default_PRIMER_CONC = 5 / 10000000 #uM

# Ionic defaults:
Default_NA  = 50   # mM
Default_MG  = 1.5  # mM
Default_dNTP = 0.8 # nM

# Multiplex Defaults:
MINIMUM_SPACING = 50
MELT_TM_TOLERANCE = 1.0
MAX_CROSS_DIMERS = 1000.0

# Input file default
DEFAULT_EXCEL_PATH = 'data/blood_borne_pathogens.xlsx' # replace with desired file path