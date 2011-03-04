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

import sys,os,shutil
from subprocess import *
import pyfits
import ReduceLog

vnum = "1.5"
#sd = '/nfs/atapplic/malt/reduce/'
sd = '/epp/atapplic/malt/malt90-analysis-code/reduce/' #Seems to be new location
data_dir = '/DATA/MALT_1/MALT90/data/'

def main():
	source = sys.argv[1]
	if source[0] == "G": #This is an individual source
		sources = [source]
	else: #Accept a date
		redlog = ReduceLog.ReduceLog()
		sources = redlog.find_files_with_date(source)
	print(sources)
	for one_source in sources:
		do_reduction(one_source)

def do_reduction(source):
	lines,freqs = setup_lines()
	try:
		override = sys.argv[2].strip()
		print(override)
	except IndexError:
		override = "none"
		print("No override command")
	over_ldata = False
	over_gzilla_ind = False
	over_gzilla_both = False
	if override == "ldata":
		over_ldata = True
	elif override == "gzilla_ind":
		over_gzilla_ind = True
	elif override == "gzilla_both":
		over_gzilla_both = True
	filenames = [source+'_GLat.rpf',source+'_GLon.rpf']
	do_livedata(filenames,lines,over_ldata)
	filenames = [source+'_GLat.sdfits',source+'_GLon.sdfits']
	do_gridzilla(source,filenames,lines,freqs,over_gzilla_ind,over_gzilla_both)
	create_source_folder(source,lines)

def setup_lines():
	### Malt90 Main Survey ###
#	lines = ["n2hp","13cs","h41a","ch3cn",\
	#	 "hc3n","13c34s","hnc","hc13ccn",\
	 #        "hcop","hcn","hnco413","hnco404",\
	  #       "c2h","hn13c","sio","h13cop"]
	#freqs = [93173.772, 92494.303, 92034.475, 91985.316, \
		# 90978.989, 90926.036, 90663.572, 90593.059, \
		 #89188.526, 88631.847, 88239.027, 87925.238, \
		 #87316.925, 87090.850, 86847.010, 86754.330]	
	### Malt90 backward sources (first day) st###
	lines = ["ch3cn","h41a","13cs","n2hp",
		 "hc13ccn","hnc","13c34s","hc3n",
		 "hnco404","hnco413","hcn","hcop",
		 "c2h","hn13c","sio","h13cop"]

	freqs = [91985.316, 92034.475,  92494.303, 93173.772,
		 90593.059, 90663.572,  90926.036, 90978.989,
		 87925.238, 88239.027,  88631.847, 89188.526,
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

def do_livedata(filenames,lines,over_ldata):
	make_dirs("livedata",lines)
	redlog = ReduceLog.ReduceLog()
	for filename in filenames:
		print(filename)
		if os.path.exists(data_dir+'renamed/'+filename):
			#print(filename)
			ftemp = filename
			ldata_done = redlog.check_val(ftemp.replace('.rpf',''),"ldata",vnum)
			#print("Ldata")
			#print(ldata_done)
			if ldata_done and not over_ldata:
				p = Popen(["glish",'-l',sd+'ldata_malt90.g','-plain',filename])
				p.wait()
				redlog.set_val(ftemp.replace('.rpf',''),"ldata",vnum)
		

def do_gridzilla(source,filenames,lines,freqs,over_gzilla_ind,over_gzilla_both):
	print(over_gzilla_both)
	make_dirs("gridzilla",lines)
	redlog = ReduceLog.ReduceLog()
	if len(filenames) == 2:
		fn1 = filenames[0]
		fn2 = filenames[1]
		gzilla_do1 = redlog.check_val(fn1.replace(".sdfits",""),"gzilla",vnum) 
		gzilla_do2 = redlog.check_val(fn2.replace(".sdfits",""),"gzilla",vnum) 
	for fn in filenames:
		fail_flag = False
		gzilla_do = redlog.check_val(fn.replace(".sdfits",""),"gzilla",vnum) 
		for i,(line,freq) in enumerate(zip(lines,freqs)):
			filein  = fn.replace(".sdfits","")+'_'+line+'.sdfits'
			fileout = fn.replace(".sdfits","")+'_'+line
			#print(filein)
			if gzilla_do or over_gzilla_ind:
				q = Popen(["glish",'-l',sd+'gzill_malt90.g',\
					line,str(freq),str(i),fileout,filein])
				q.wait() 
				try:
					shutil.move(data_dir+'/gridzilla/'+line+'/'+fileout+'.beamRSS.fits',data_dir+'/gridzilla/'+line+'/'+'flotsam')
					shutil.move(data_dir+'/gridzilla/'+line+'/'+fileout+'.spectracounts.fits',data_dir+'/gridzilla/'+line+'/'+'flotsam')
				except IOError:
					fail_flag = True
		if not fail_flag:
			redlog.set_val(fn.replace(".sdfits",""),"gzilla",vnum)
	if len(filenames) == 2:
		for i,(line,freq) in enumerate(zip(lines,freqs)):
			file1   = fn1.replace(".sdfits","")+'_'+line+'.sdfits'
			file2   = fn2.replace(".sdfits","")+'_'+line+'.sdfits'
			fileout = source+'_'+line
			if gzilla_do1 or gzilla_do2 or over_gzilla_both:
				q = Popen(["glish",'-l',sd+'gzill_malt90.g',\
					line,str(freq),str(i),fileout,file1,file2])
				q.wait()
				try:
					shutil.move(data_dir+'/gridzilla/'+line+'/'+fileout+'.beamRSS.fits',data_dir+'/gridzilla/'+line+'/'+'flotsam')
					shutil.move(data_dir+'/gridzilla/'+line+'/'+fileout+'.spectracounts.fits',data_dir+'/gridzilla/'+line+'/'+'flotsam')
				except IOError:
					print("Failed to make joint source")
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
	
if __name__ == '__main__':
	main()

