#This file stores variables for use 
#throughout the MALT90 reduction pipeline

#moment_map
#base_data = "/DATA/MALT_1/MALT90/data/"
#Redundant with data_dir

base = '/DATA/MALT_1/MALT90/'

#reduce_malt
#vnum = "1.6"
vnum = {"rename":"1.5","ldata":"1.6","gzilla":"1.6","arrange":"1.6","mommaps":"1.6"}
#sd = '/nfs/atapplic/malt/reduce/'                                                                                                                                                        
sd = '/epp/atapplic/malt/malt90-analysis-code/reduce/' #Seems to be new location
data_dir = base+'data/'
data_dir_y1 = base+'data/'+'year1/'
data_dir_y2 = base+'data/'+'year2/'

#ReduceLog
log_location = base+'reduction_log.txt'
lock    = base+'lock.txt'
source_dir = base+'raw_data/' #This is rather badly named
#classification_file = sd+'Malt90Catalog_classifications_v16.txt'
#classification_file = sd+'MALT90_classifications_2011.txt'
classification_file = sd+'MALT90_classifications_2012.txt'

#preprocess
rename_dir = '/DATA/MALT_1/MALT90/data/renamed/'
#data_dir = '/DATA/MALT_1/MALT90/raw_data/'

#examine_cal

#data_dir = "/DATA/MALT_1/MALT90/cal/"
cal_dir = base+"cal/"
ver_dir = base+"data/verification/"

#These are paths used by do_by_hand.py
byhand_rename_dir = base+'data/byhand/renamed/'
byhand_data_dir   = base+'data/byhand/'
