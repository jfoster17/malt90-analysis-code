#!/usr/bin/env python
# encoding: utf-8
"""
preprocess_malt.py

Summary:
Step one of the Malt90 reduction pipeline
-Renames files based on source
-Tags bad (too small) files in log with B on end of sources name

Commands:
preprocess_malt.py 2010-07-29_0646-M516.rpf
Do a single file, ignorning/overwriting if it has been done before.

preprocess_malt.py 2010-07-29 
Do all unprocessed files from a given date

preprocess_maly.py all
Do all unrpocessed files

A second argument of "obs" will denote that observations are currently
underway. Too-small files will be ignored so they are not falsely
entered as bad.

Changes:
Version 1.4 disables ASAP smoothing because it was corrupting the velocity axis. 
Also changes directories to preserve old reduction for comparison
Version 1.5 correctly identifies SiO and HN13C lines data (they were swapped in previous versions)

"""
import sys, os, glob, shutil
from asap import *
import ReduceLog

sd = '/nfs/atapplic/malt/reduce/'
rename_dir = '/DATA/MALT_1/MALT90/data/renamed/'
data_dir = '/DATA/MALT_1/MALT90/raw_data/'
vnum = "1.5"

def main():
	argument = None
	in_middle_of_obs = None
	try: 
		argument = sys.argv[1]
	except:
		pass
	try:
		in_middle_of_obs = sys.argv[2]
	except:
		pass
	if in_middle_of_obs:
		in_middle_of_obs = True
	if argument.endswith('all'): #Get all undone files
		files_to_process = get_new_files(in_middle_of_obs=in_middle_of_obs)
	elif argument.endswith('.rpf'): #A single file
		files_to_process = [argument]
	else: #Assume it is a date
		files_to_process = get_new_files(argument,in_middle_of_obs = in_middle_of_obs)
	rename_files(files_to_process)

def get_new_files(date = "all",in_middle_of_obs=False):
	"""Get list of files to update"""
	redlog = ReduceLog.ReduceLog()
	redlog.update(date,in_middle_of_obs=in_middle_of_obs)
	if date != "all":
		datestring = date
	else:
		datestring = "20" #This will always be in a date
	return(redlog.find_new_files(vnum,datestring))	

def rename_files(filelist):
	"""Load files into ASAP. Smooth references. Lookup name/check size and rename"""
	redlog = ReduceLog.ReduceLog()
	for new_file in filelist:
		s = scantable(data_dir+new_file,average=False)	
		id,source_for_log = redlog.get_name(new_file)
		renamed_file = source_for_log+".rpf"
		print("Saving "+new_file+" as "+source_for_log+"...")
		#s.save(smoothdir+source_for_log+'.sdfits','SDFITS',overwrite=True)
		shutil.copyfile(data_dir+new_file,rename_dir+source_for_log+".rpf")
		
		redlog.set_val(source_for_log,"rename",vnum)
		print("")


if __name__ == '__main__':
	main()


