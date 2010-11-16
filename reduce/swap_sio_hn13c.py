#!/usr/bin/env python
# encoding: utf-8
"""
swap_sio_hn13c.py
"""

import sys
import os
import pyfits

data_dir = "/mako3/MALT_1/MALT90/data/"

def main():
	"""HN13C and SiO were swapped in the v1.4 data pipeline (and prior). This script renames files to fix this."""
	
	fix_gridzilla()
	#fix_moment_maps()
	#fix_sources()
	
	pass

def fix_headers(subfolder,currIF,newIF):
	"""Update headers to reflect hn13c/sio swap"""
	all_files = os.listdir(data_dir+subfolder)
	#all_files = ["test.fits"]
	for curr_file in all_files:
		hdulist = pyfits.open(data_dir+subfolder+curr_file,mode='update')
		prihdr = hdulist[0].header
		for i,entry in enumerate(prihdr.ascardlist()):
			blah = str(prihdr[i]).find(currIF)
			if blah != -1:
				prihdr[i] = str(prihdr[i]).replace(currIF,newIF)
		prihdr.add_history("Fixed header due to hn13c/sio swap")
		hdulist.close()

def fix_gridzilla():
	"""Fix the gridzilla folder, including headers"""
	fix_headers("gridzilla/hn13c/","hn13c","sio")
	all_files = os.listdir(data_dir+"gridzilla/hn13c/")
	for curr_file in all_files:
		new_file = curr_file.replace("hn13c","sio")
		os.rename(data_dir+"gridzilla/hn13c/"+curr_file,data_dir+"gridzilla/hn13c/"+new_file)
		
	fix_headers("gridzilla/sio/","sio","hn13c")
	all_files = os.listdir(data_dir+"gridzilla/sio/")
	for curr_file in all_files:
		new_file = curr_file.replace("sio","hn13c")
		os.rename(data_dir+"gridzilla/sio/"+curr_file,data_dir+"gridzilla/sio/"+new_file)
		
	#This is the basic switcheroo
	os.rename(data_dir+"gridzilla/hn13c/",data_dir+"gridzilla/temp")
	os.rename(data_dir+"gridzilla/sio/",data_dir+"gridzilla/hn13c")
	os.rename(data_dir+"gridzilla/temp/",data_dir+"gridzilla/sio")
	
		
#def fix_moment_maps():
#	"""Fix moment_maps folder, including headers"""
#	all_dirs = os.listdir(data_dir+"moment_maps/hn13c/")
#	for curr_directory in all_dirs:
#		all_files = os.listdir(data_dir+"moment_maps/hn13c/"+curr_directory)
		
#def fix_sources():
#	"""Fix links in sources and moment_maps as in fix_moment_maps"""

if __name__ == '__main__':
	main()

