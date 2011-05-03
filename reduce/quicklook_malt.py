#!/usr/bin/env python
# encoding: utf-8
"""
quicklook_malt.py

A quick reduction of hco+ (and/or hnc) to check that we actually
see a source.

Call with the name of your source and which direction you want
to look at (leave off direction to combine both)
quicklook_malt -s G320.020+00.100 -d GLat

That is, source name should follow -s and
direction (GLat or GLon) should follow -d

"""

import sys,os,getopt
import pyfits
import numpy as np
import preprocess_malt,reduce_malt,ReduceLog
import datetime

def main():
	try:
		opts,args = getopt.getopt(sys.argv[1:],"s:d:")
	except getopt.GetoptError,err:
		print str(err)
		print __doc__
		sys.exit(2)
	direction = None
	for o,a in opts:
		if o == "-s":
			source = a
		if o == "-d":
			direction = a

	#Make sure that log is up-to-date.
	#This step is quick, so just do +/-
	#on current UT. 
	now_utc = datetime.datetime.utcnow()
	today = datetime.date(now_utc.year,now_utc.month,now_utc.day)
	plus1 = today + datetime.timedelta(days=1)
	minu1 = today - datetime.timedelta(days=1)

	for dates in [today,plus1,minu1]:
		files_to_process = preprocess_malt.get_new_files(dates.isoformat(),in_middle_of_obs = True)
		preprocess_malt.rename_files(files_to_process)
		
	#redlog = ReudceLog.ReduceLog()

	reduce_malt.do_reduction(source,force_list=['ldata','gzilla','mommaps'],quicklook=True,onlyone=direction)

if __name__ == '__main__':
	main()

