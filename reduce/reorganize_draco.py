#!/usr/bin/env python
# encoding: utf-8
"""
reorganize_mako.py

A script to reorganize Malt90 data on Draco.
"""

import sys,os,shutil
from subprocess import *
import pyfits

vnum = "1.5"
data_dir = '/DATA/MALT_1/MALT90/data/'

def main():
	all_files = os.listdir(data_dir+"gridzilla/13c34s/")
	sources = []
	for filename in all_files:
	#	print(filename)
		if (filename.find("GLat") < 0) and (filename.find("GLon") < 0 ) and filename.startswith("G"):
			print(filename)
			sources.append(filename.partition("_")[0])
	print(sources)
	for one_source in sources:
			do_reduction(one_source)

def do_reduction(source):
	lines,freqs = setup_lines()
	create_source_folder(source,lines)

def setup_lines():
	### Malt90 Main Survey ###
	lines = ["n2hp","13cs","h41a","ch3cn",\
		 "hc3n","13c34s","hnc","hc13ccn",\
	         "hcop","hcn","hnco413","hnco404",\
	         "c2h","hn13c","sio","h13cop"]
	freqs = [93173.772, 92494.303, 92034.475, 91985.316, \
		 90978.989, 90926.036, 90663.572, 90593.059, \
		 89188.526, 88631.847, 88239.027, 87925.238, \
		 87316.925, 87090.850, 86847.010, 86754.330]	
	return(lines,freqs)

def make_dirs(dirname,lines):
	try:
		os.mkdir(data_dir+dirname)
	except OSError:
		pass
	for line in lines:
		try:
			os.mkdir(data_dir+dirname+"/"+line)
		except OSError:
			pass

def create_source_folder(source,lines):
	"""Create a folder of symlinks for each source"""
	try:
		os.mkdir(data_dir+"sources")
	except OSError:
		pass
	for line in lines:
		try:
			os.mkdir(data_dir+'sources/'+source)
		except OSError:
			pass
		filename = source+'_'+line+'_MEAN.fits'
		src = data_dir+'gridzilla/'+line+'/'+filename
		targ = data_dir+'sources/'+source+'/'+filename
		try:
			os.unlink(targ)
		except OSError:
			pass
		try:
			os.symlink(src,targ)
		except OSError:
			pass
		momsrc = data_dir+'moment_maps/'+line+'/'+source+'_'+line+'_MEAN_mom_maps'
		momtarg = data_dir+'sources/'+source+'/'+source+'_'+line+'_MEAN_mom_maps'
		#Uncomment to restore moment map copying
		try:
			shutil.rmtree(momtarg)
		except OSError:
			pass
		try:
			shutil.copytree(momsrc,momtarg)
		except OSError:
			pass
	
if __name__ == '__main__':
	main()

