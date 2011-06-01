#!/usr/bin/env python
# encoding: utf-8
"""
Copy over all the year 1 stuff
"""

import sys
import os
import glob
import reduce_malt
import malt_params
import shutil

def main():
	f = open('malt90_velocities_year1.txt','r')
	source_list = []
	for line in f:
		name,vhand = line.split()
		source_list.append(name)
	try:
		reverse = sys.argv[1]
	except:
		reverse = False
	#source_list = ['G318.725-00.224'] #Ugly code to do one source
	for sourcename in source_list:
		copy_year1_source("gridzilla",sourcename,reverse)
		copy_year1_source("livedata",sourcename,reverse)
		copy_year1_source("mommaps",sourcename,reverse)
		copy_year1_source("mommaps",sourcename,reverse)
		copy_year1_source_plain("renamed",sourcename,reverse)
		copy_year1_source_plain("sources",sourcename,reverse)

def copy_year1_source(directory,sourcename,reverse):
	lines,freqs,ifs = reduce_malt.setup_lines()
	for line in lines:
		if reverse:
			search_dir = os.path.join(malt_params.data_dir,
						  'year1',directory,line)
		else:
			search_dir = os.path.join(malt_params.data_dir,directory,line)
		print(search_dir+'/'+sourcename+'*')
		files = glob.glob(search_dir+'/'+sourcename+'*')
		for fin in files:
			if reverse:
				fout = fin.replace('/data/year1/','/data/')
			else:
				fout = fin.replace('/data/','/data/year1/')
			print("From: "+fin)
			print("To: "+fout)
	     		shutil.move(fin,fout)

def copy_year1_source_plain(directory,sourcename,reverse):
	if reverse:
		search_dir = os.path.join(malt_params.data_dir,
					  'year1',directory)
	else:
		search_dir = os.path.join(malt_params.data_dir,directory)
	print(search_dir+'/'+sourcename+'*')
	files = glob.glob(search_dir+'/'+sourcename+'*')
	for fin in files:
		if reverse:
			fout = fin.replace('/data/year1/','/data/')
		else:
			fout = fin.replace('/data/','/data/year1/')
		print("From: "+fin)
		print("To: "+fout)
		shutil.move(fin,fout)


if __name__ == '__main__':
	main()

