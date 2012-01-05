#!/usr/bin/env python
# encoding: utf-8
"""
Copy over all the year 2 stuff
"""

import sys
import os
import glob
import reduce_malt
import malt_params
import shutil

dest = 'year2'


def main():
	f = open('malt90_velocities_year2.txt','r')
	source_list = []
	for line in f:
		name,vhand = line.split()
		source_list.append(name)
	try:
		reverse = sys.argv[1]
	except:
		reverse = False
	#source_list = ['G007.276-00.531']
	#source_list = ['G318.725-00.224'] #Ugly code to do one source
	for sourcename in source_list:
		copy_source("gridzilla",sourcename,reverse)
		copy_source("livedata",sourcename,reverse)
		copy_source("mommaps",sourcename,reverse)
		copy_source_plain("renamed",sourcename,reverse,False)
		copy_source_plain("sources",sourcename,reverse,True)

def copy_source(directory,sourcename,reverse):
	lines,freqs,ifs = reduce_malt.setup_lines()
	for line in lines:
		if reverse:
			search_dir = os.path.join(malt_params.data_dir,
						  dest,directory,line)
		else:
			search_dir = os.path.join(malt_params.data_dir,directory,line)
		#print(search_dir+'/'+sourcename+'*')
		files = glob.glob(search_dir+'/'+sourcename+'*')
		for fin in files:
			if reverse:
				fout = fin.replace('/data/'+dest+'/','/data/')
			else:
				fout = fin.replace('/data/','/data/'+dest+'/')
			#print("From: "+fin)
			#print("To: "+fout)
	     		#shutil.move(fin,fout)

def copy_source_plain(directory,sourcename,reverse,symlink):
	if reverse:
		search_dir = os.path.join(malt_params.data_dir,
					  dest,directory)
	else:
		search_dir = os.path.join(malt_params.data_dir,directory)
	print(search_dir+'/'+sourcename+'*')
	files = glob.glob(search_dir+'/'+sourcename+'*')
	for fin in files:
		if reverse:
			fout = fin.replace('/data/'+dest+'/','/data/')
			#fout = fout.replace('/DATA/MALT_1/MALT90/data/','../../')
		else:
			fout = fin.replace('/data/','/data/'+dest+'/')
			#if symlink:
			#	fout = fout.replace('/DATA/MALT_1/MALT90/data/','../../')
		print("From: "+fin)
		print("To: "+fout)
		shutil.move(fin,fout)
	if symlink:
		print(fout+'/*MEAN*')
		files = glob.glob(fout+'/*MEAN*')
		print(files)
		for fsin in files:
			fsout = fsin.replace('/DATA/MALT_1/MALT90/data/year2/sources/','../../gridzilla/')
			print("From symlink:"+fsin)
			print("To   symlink:"+fsout)
			#shutil.move(fsin,fsout)
			try:
				os.unlink(fsin)
			except OSError:
				pass
			try:
				os.symlink(fsout,fsin)
			except OSError:
				pass

if __name__ == '__main__':
	main()

