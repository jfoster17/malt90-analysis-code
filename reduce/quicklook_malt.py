#!/usr/bin/env python
# encoding: utf-8
"""
quicklook_malt.py is a quick reduction of hco+ and hnc to 
check that we actually see a source.

Call with the name of your source and which direction you want
to exaimne (leave off direction to combine both GLat and GLon)
quicklook_malt.py -s G320.020+00.100 -d GLat

Options
-s : Source The full name of the source you would like to reduce.
-d : Direction The scan direction (GLat or GLon) to reduce. 
     Leave blank to reduce both (reducing both maps takes longer).
-h : Display this help

Verification images are created in /DATA/MALT_1/MALT90/data/verification/
"""

import sys,os,getopt
import pyfits
import numpy as np
import datetime
import malt_params as malt
import idl_stats
import preprocess_malt,reduce_malt,ReduceLog

def main():
	try:
		opts,args = getopt.getopt(sys.argv[1:],"s:d:h")
	except getopt.GetoptError,err:
		print(str(err))
		print(__doc__)
		sys.exit(2)
	direction = None
	for o,a in opts:
		if o == "-s":
			source = a
		if o == "-d":
			direction = a
		if o == "-h":
			print(__doc__)
			sys.exit(1)
	#Make sure that log is up-to-date.
	#This step is quick, so just do +/- on current UT. 
	now_utc = datetime.datetime.utcnow()
	today = datetime.date(now_utc.year,now_utc.month,now_utc.day)
	plus1 = today + datetime.timedelta(days=1)
	minu1 = today - datetime.timedelta(days=1)

	for dates in [today,plus1,minu1]:
		files_to_process = preprocess_malt.get_new_files(
			           dates.isoformat(),in_middle_of_obs = True)
		preprocess_malt.rename_files(files_to_process)
		
	reduce_malt.do_reduction(source,
				 force_list=['ldata','gzilla','mommaps'],
				 quicklook=True,onlyone=direction)
#	print("This is original direction")
	#print(direction)
	#reduce_malt.make_verification_plots(source,direction=direction)

		

if __name__ == '__main__':
	main()

