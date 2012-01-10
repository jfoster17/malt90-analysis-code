#!/usr/bin/env python
# encoding: utf-8
"""
General script to copy a given list os sources from/to specified folders.
-l : file from which to pull list of sources to move
-o : origin folder. To move from pipeline use 'data'. From year2 use 'data/year2'
-d : destination folder; same format as origin
-t : test run only (print files to move)
-h : display this help

For example, to move the appropriate files from the reduction pipeline to the year2
folder I did
python copy_sources.py -l year2_plain_sources.txt -d data/year2 -o data
python copy_sources.py -l year2_byhand_sources.txt -d data/year2 -o data/byhand

I did run this with the livedata and renamed folders commented out. But for 
long-term storage of files I think they should be uncommented and the above 
commands run again.

Also, now that I have a good procedure for by-hand data I should really do
year1 byhand and make better byhand notes.
"""

import sys,os,glob,getopt
import reduce_malt
import malt_params
import shutil

def main():
	try:
		opts,args=getopt.getopt(sys.argv[1:],"l:o:d:th")
	except getopt.GetoptError,err:
		print(str(err))
		print(__doc__)
		sys.exit(2)
	test_run = False
	for o,a in opts:
		if o == "-l":
			source_list_file = a
		if o == "-o":
			origin = a
		if o == "-d":
			destination = a
		if o == "-t":
			test_run = True
		if o == "-h":
			print(__doc__)
			sys.exit(1)
	f = open(source_list_file,'r')
	source_list = []
	for line in f:
		name = line.strip()
		source_list.append(name)
	for sourcename in source_list:
		copy_source("gridzilla",sourcename,origin,destination,test_run)
#		copy_source("livedata",sourcename,origin,destination,test_run)
		copy_source("mommaps",sourcename,origin,destination,test_run)
	#	copy_source_plain("renamed",sourcename,origin,destination,test_run,False)
		copy_source_plain("sources",sourcename,origin,destination,test_run,True)
		clean_links(sourcename,destination)

def copy_source(directory,sourcename,origin,destination,test_run):
	lines,freqs,ifs = reduce_malt.setup_lines()
	for line in lines:
		search_dir = os.path.join(malt_params.base,origin,directory,line)
#		print("Search Dir=")
	#	print(search_dir)
		files = glob.glob(search_dir+'/'+sourcename+'*')
		for fin in files:
			if ((("_2" not in fin) and ("_3" not in fin)) or ("_2" in sourcename) or ("_3" in sourcename)): #Logic to remove unwanted _2
				fout = fin.replace(origin,destination)
				print("From: "+fin)
				print("To: "+fout)
				if not test_run:
					try:
						shutil.rmtree(fout)
					except OSError:
						pass
					try:
						os.remove(fout)
					except OSError:
						pass
					shutil.move(fin,fout)

def clean_links(sourcename,destination):
	search_dir = os.path.join(malt_params.base,
				  destination,"sources")
	files = glob.glob(search_dir+'/'+sourcename+'/'+'*MEAN.fits')
	for fsin in files:
	#	print(fsin)
		line = os.path.basename(fsin).split("_")[1]
		fsout = "../../gridzilla/"+line+'/'+os.path.basename(fsin)
#		print(newline)
#		fsout = fsin.replace(search_dir,'../../gridzilla/')
		try:
			os.unlink(fsin)
		except OSError:
			pass
		try:
			os.symlink(fsout,fsin)
		except OSError:
			pass

def copy_source_plain(directory,sourcename,origin,destination,test_run,symlink):
	search_dir = os.path.join(malt_params.base,
					  origin,directory)
	
	files = glob.glob(search_dir+'/'+sourcename+'*')
	for fin in files:
		if ((("_2" not in fin) and ("_3" not in fin)) or ("_2" in sourcename) or ("_3" in sourcename)): #Logic to remove unwanted _2
			fout = fin.replace(origin,destination)
			print("From: "+fin)
			print("To: "+fout)
			if not test_run:
				try:
					shutil.rmtree(fout)
				except OSError:
					pass
				try:
					os.remove(fout)
				except OSError:
					pass
				shutil.move(fin,fout)

if __name__ == '__main__':
	main()

