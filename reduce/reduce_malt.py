#!/usr/bin/env python
# encoding: utf-8
"""
reduce_malt.py

A script to reduce Malt90 data.
Executes glish scripts which control livedata/gridzilal
Raw-Data->renamed/source
Livedata->livedata/line
Gridzilla->gridzilla/line
Version 1.5 correctly identifies SiO and HN13C lines data (they were swapped in previous versions)
"""

import sys,os,shutil,getops
from subprocess import *
import pyfits
import ReduceLog

vnum = "1.5"
#sd = '/nfs/atapplic/malt/reduce/'
sd = '/epp/atapplic/malt/malt90-analysis-code/reduce/' #Seems to be new location
data_dir = '/DATA/MALT_1/MALT90/data/'

def main():
	try:
		opts,args=getopt.getopt(sys.argv[1:], "nsafi")
	except getopt.GetoptError,err:
		print str(err)
		usage()
		sys.exit(2)
	force_list  = []
	ignore_list = []
	do_source = False
	do_date   = False
	do_all    = False

	for o,a in opts:
		if o == "-f":
			force_list = list(a)
		elif o == "-i":
			ignore_list = list(a)
		elif o == "-n":
			do_date = True
			date = a
		elif o == "-s":
			do_source = True
			source = a
		elif o == "-a":
			do_all = True
		else:
			assert False, "unhandled option"
	redlog = ReduceLog.ReduceLog()
	if do_date:
		sources = redlog.find_files_with_date(source)
	elif do_source:
		sources = [source]
	elif do_all:
		sources = redlog.find_undone(source,vnum)
	print(sources)
	for one_source in sources:
		do_reduction(one_source,force_list,ignore_list)

def do_reduction(source,force_list=None,ignore_list=None):
	lines,freqs = setup_lines()
	### Do Livedata ###
	filenames = [source+'_GLat.rpf',source+'_GLon.rpf']
	if 'ldata' in force_list:
		do_livedata(filenames,lines,force=True)
	elif 'ldata' in ignore_list:
		pass
	else:
		do_livedata(filenames,lines,force=False)
	### Do Gridzilla ###
	filenames = [source+'_GLat.sdfits',source+'_GLon.sdfits']
	if 'gzilla' in force_list:
		do_gridzilla(source,filenames,lines,freqs,force=True)
	elif 'gzilla' in ignore_list:
		pass
	else:
		do_gridzilla(source,filenames,lines,freqs,force=False)
	### Do Reorganization ###
	### I think it is fine to always do this step ###
	create_source_folder(source,lines)

	### Do moment maps ###
	if 'mommaps' in force_list:
		do_mommaps(source,filenames,force=True)
	elif 'mommaps' in ignore_list:
		pass
	else:
		do_mommaps(source,filenames,force=False)

def setup_lines():
	### Malt90 Main Survey ###
	lines = ["n2hp","13cs","h41a","ch3cn",
		 "hc3n","13c34s","hnc","hc13ccn",
	         "hcop","hcn","hnco413","hnco404",
	         "c2h","hn13c","sio","h13cop"]
	freqs = [93173.772, 92494.303, 92034.475, 91985.316, 
		 90978.989, 90926.036, 90663.572, 90593.059, 
		 89188.526, 88631.847, 88239.027, 87925.238, 
		 87316.925, 87090.850, 86847.010, 86754.330]	
	### Malt90 backward sources (first day) st###
#	lines = ["ch3cn","h41a","13cs","n2hp",
#		 "hc13ccn","hnc","13c34s","hc3n",
#		 "hnco404","hnco413","hcn","hcop",
#		 "c2h","hn13c","sio","h13cop"]

#	freqs = [91985.316, 92034.475,  92494.303, 93173.772,
#		 90593.059, 90663.572,  90926.036, 90978.989,
#		 87925.238, 88239.027,  88631.847, 89188.526,
#		 87316.925, 87090.850, 86847.010, 86754.330]	

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

def do_livedata(filenames,lines,force=False):
	make_dirs("livedata",lines)
	redlog = ReduceLog.ReduceLog()
	for filename in filenames:
		print(filename)
		if os.path.exists(data_dir+'renamed/'+filename):
			#print(filename)
			ftemp = filename
			ldata_needed = redlog.check_val(ftemp.replace('.rpf',''),"ldata",vnum)
			#print("Ldata")
			#print(ldata_done)
			if ldata_needed or force:
				p = Popen(["glish",'-l',sd+'ldata_malt90.g','-plain',filename])
				p.wait()
				redlog.set_val(ftemp.replace('.rpf',''),"ldata",vnum)
		

def do_gridzilla(source,filenames,lines,freqs,force=False):
	# print(over_gzilla_both)
	make_dirs("gridzilla",lines)
	redlog = ReduceLog.ReduceLog()
	if len(filenames) == 2:
		fn1 = filenames[0]
		fn2 = filenames[1]
		gzilla_do1 = redlog.check_val(fn1.replace(".sdfits",""),"gzilla",vnum) 
		gzilla_do2 = redlog.check_val(fn2.replace(".sdfits",""),"gzilla",vnum) 
	for fn in filenames:
		fail_flag = False
		gzilla_needed = redlog.check_val(fn.replace(".sdfits",""),"gzilla",vnum) 
		for i,(line,freq) in enumerate(zip(lines,freqs)):
			filein  = fn.replace(".sdfits","")+'_'+line+'.sdfits'
			fileout = fn.replace(".sdfits","")+'_'+line
			#print(filein)
			if gzilla_needed or force:
				q = Popen(["glish",'-l',sd+'gzill_malt90.g',\
					line,str(freq),str(i),fileout,filein])
				q.wait() 
				try:
					shutil.move(data_dir+'/gridzilla/'+line+'/'+fileout+'.beamRSS.fits',data_dir+'/gridzilla/'+line+'/'+'flotsam/')
					shutil.move(data_dir+'/gridzilla/'+line+'/'+fileout+'.spectracounts.fits',data_dir+'/gridzilla/'+line+'/'+'flotsam/')
				except IOError:
					fail_flag = True
					print("Failed to move flotsam files")
		if not fail_flag:
			redlog.set_val(fn.replace(".sdfits",""),"gzilla",vnum)
	if len(filenames) == 2:
		for i,(line,freq) in enumerate(zip(lines,freqs)):
			file1   = fn1.replace(".sdfits","")+'_'+line+'.sdfits'
			file2   = fn2.replace(".sdfits","")+'_'+line+'.sdfits'
			fileout = source+'_'+line
			if gzilla_do1 or gzilla_do2 or force:
				q = Popen(["glish",'-l',sd+'gzill_malt90.g',\
					line,str(freq),str(i),fileout,file1,file2])
				q.wait()
				try:
					shutil.move(data_dir+'/gridzilla/'+line+'/'+fileout+'.beamRSS.fits',data_dir+'/gridzilla/'+line+'/'+'flotsam/')
					shutil.move(data_dir+'/gridzilla/'+line+'/'+fileout+'.spectracounts.fits',data_dir+'/gridzilla/'+line+'/'+'flotsam/')
				except IOError:
					print("Failed to make joint source")
					#This does not return a useful error
def create_source_folder(source,lines):
	"""Create a folder of symlinks for each source"""
	redlog = ReduceLog.ReduceLog()
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
			os.symlink(src,targ)
		except OSError:
			pass
		momsrc = data_dir+'moment_maps/'+line+'/'+source+'_'+line+'_MEAN_mom_maps'
		momtarg = data_dir+'sources/'+source+'/'+source+'_'+line+'_MEAN_mom_maps'
		try:
			shutil.copytree(momsrc,momtarg)
		except OSError:
			pass
		try:
			hdulist = pyfits.open(targ,mode='update')
			prihdr = hdulist[0].header
			prihdr.update('M90PIPEV',vnum,'Malt90 Pipeline Version')
			prihdr.update('BUNIT','K','Antenna Temperature')
			hdulist.flush()
		except:
			print("Header update failed")
	files_involved = [source+"_GLat",source+"_GLon"]
	for file_involved in files_involved:
		#print(file_involved)
		try:
			redlog.set_val(file_involved,"arrange",vnum)
		except:
			print("Failed to update log for "+file_involved+". Maybe it does not exist?")
	
def do_mommaps(source,lines):
	"""Make moment maps for source and place them correctly
	Determine a source velocity (currently from table)
	Make moment map for GLon, GLat, combined in mommaps/line
	Copy the combined moment map into sources/ folder
	"""

	redlog = ReduceLog.ReduceLog()
	files_involved = [source+"_GLat",source+"_GLon"]
	print("I would be doing a moment map")
	for file_involved in files_involved:
		redlog.set_val(file_involved,"mommaps",vnum)
	pass

if __name__ == '__main__':
	main()

