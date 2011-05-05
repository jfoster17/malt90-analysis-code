#!/usr/bin/env python
# encoding: utf-8
"""
reduce_malt.py is the main script to reduce Malt90 data.
It executes glish scripts which control custom versions of
livedata/gridzila and places data in the following locations:
Raw-Data->renamed/source
Livedata->livedata/line
Gridzilla->gridzilla/line

Typically this script is not called directly, and is only
called by quicklook_malt.py and do_night_maly.py. The command
line options exist mainly to facilitate re-reduction of 
sources which failed to reduce the first time.

Options
-n : Night  -- reduce all files taken on a give data (YYYY-MM-DD)
-s : Source -- reduce a given source (GXXX.XXX-XX.XXX)
-a : All    -- reduce all files reduced under old versions of the pipeline
               this should, obviously, be used sparingly. In particular
               version 1.4 and 1.5 are functionarly identical, since an
	       external script was used to swap the SiO and HNC13C lines.
-f : Force  -- force the reduction of specific steps even if marked as done
               in the reduction log. Enter as a comma-separated list.
	       Valid options: ldata,gzilla,mommaps
-i : Ignore -- ignore the listed steps when reducing, even if the version
               listed in the reduction log is out-of-date. Will crash
	       if previous steps do not exists. Enter as comma-separated list
	       Valid options: ldata,gzilla,mommaps

--- Changelog ---
Version 1.5 correctly identifies SiO and HN13C lines data (they were swapped in previous versions)
"""

import sys,os,shutil,getopt
from subprocess import *
import pyfits
import ReduceLog
import moment_map
import malt_params as malt

def main():
	try:
		opts,args=getopt.getopt(sys.argv[1:], "n:s:af:i:")
	except getopt.GetoptError,err:
		print str(err)
		print(__doc__)
		sys.exit(2)
	force_list  = []
	ignore_list = []
	do_source = False
	do_date   = False
	do_all    = False

	for o,a in opts:
		if o == "-f":
			force_list = a.split(',')
			print(force_list)
		elif o == "-i":
			ignore_list = a.split(',')
			print(ignore_list)
		elif o == "-n":
			do_date = True
			date = a
		elif o == "-s":
			do_source = True
			source = a
			print(source)
		elif o == "-a":
			#I think force does not work well with do_all
			do_all = True
		else:
			assert False, "unhandled option"
	redlog = ReduceLog.ReduceLog()
	if do_date:
		sources = redlog.find_files_with_date(date)
	elif do_source:
		sources = [source]
	elif do_all:
		sources = redlog.find_undone(malt.vnum)
		#print(sources)
		#print(nothing)
	try:
		print("Reducing the following source:"+sources)
	except UnboundLocalError:
		print("@@@ No sources to reduce @@@")
		print(__doc__)
		sys.exit(2)
	for one_source in sources:
		do_reduction(one_source,force_list,ignore_list)

def do_reduction(source,force_list=[],ignore_list=[],quicklook=False,onlyone=None):
	lines,freqs,ifs = setup_lines(quicklook)

	### Do Livedata ###
	if onlyone:
		filenames = [source+'_'+onlyone+'.rpf']
	else:
		filenames = [source+'_GLat.rpf',source+'_GLon.rpf']
	if 'ldata' in force_list:
		do_livedata(filenames,lines,force=True,quicklook=quicklook)
	elif 'ldata' in ignore_list:
		pass
	else:
		do_livedata(filenames,lines,force=False,quicklook=quicklook)
	### Do Gridzilla ###
	if onlyone:
		filenames = [source+'_'+onlyone+'.sdfits']
	else:
		filenames = [source+'_GLat.sdfits',source+'_GLon.sdfits']
	if 'gzilla' in force_list:
		do_gridzilla(source,filenames,lines,freqs,ifs,force=True,quicklook=quicklook)
	elif 'gzilla' in ignore_list:
		pass
	else:
		do_gridzilla(source,filenames,lines,freqs,ifs,force=False,quicklook=quicklook)
	### Do Reorganization ###
	### I think it is fine to always do this step ###
	create_source_folder(source,lines,quicklook)

	### Do moment maps ###
	if onlyone:
		filenames = [source+'_'+onlyone]
	else:
		filenames = [source+"_GLat",source+"_GLon",source]
	if 'mommaps' in force_list:
		do_mommaps(source,filenames,lines,force=True,quicklook=quicklook)
	elif 'mommaps' in ignore_list:
		pass
	else:
		do_mommaps(source,filenames,lines,force=False,quicklook=quicklook)

def setup_lines(quicklook=False):
	if quicklook:
		### Lines to use for quick look ###
		lines = ["hcop"]
		freqs = [89188.526]
		ifs   = [9]
	else:
        	### Malt90 Main Survey ###
		lines = ["n2hp","13cs","h41a","ch3cn",
			 "hc3n","13c34s","hnc","hc13ccn",
			 "hcop","hcn","hnco413","hnco404",
			 "c2h","hn13c","sio","h13cop"]
		freqs = [93173.772, 92494.303, 92034.475, 91985.316, 
			 90978.989, 90926.036, 90663.572, 90593.059, 
			 89188.526, 88631.847, 88239.027, 87925.238, 
			 87316.925, 87090.850, 86847.010, 86754.330]	
		ifs   = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16]
	### Malt90 backward sources (first day) st###
#	lines = ["ch3cn","h41a","13cs","n2hp",
#		 "hc13ccn","hnc","13c34s","hc3n",
#		 "hnco404","hnco413","hcn","hcop",
#		 "c2h","hn13c","sio","h13cop"]

#	freqs = [91985.316, 92034.475,  92494.303, 93173.772,
#		 90593.059, 90663.572,  90926.036, 90978.989,
#		 87925.238, 88239.027,  88631.847, 89188.526,
#		 87316.925, 87090.850, 86847.010, 86754.330]	

	return(lines,freqs,ifs)

def make_dirs(dirname,lines):
	try:
		os.mkdir(malt.data_dir+dirname)
	except OSError:
		pass
	for line in lines:
		try:
			os.mkdir(malt.data_dir+dirname+"/"+line)
		except OSError:
			pass

def do_livedata(filenames,lines,force=False,quicklook=False):
	make_dirs("livedata",lines)
	redlog = ReduceLog.ReduceLog()
	for filename in filenames:
		print(filename)
		if os.path.exists(malt.data_dir+'renamed/'+filename):
			#print(filename)
			ftemp = filename
			ldata_needed = redlog.check_val(ftemp.replace('.rpf',''),"ldata",malt.vnum)
			#print("Ldata")
			#print(ldata_done)
			if ldata_needed or force:
				#I think I could use check_call here to see if this dies
				if not quicklook:
					p = Popen(["glish",'-l',malt.sd+'ldata_malt90.g','-plain',filename])
					p.wait()
				else:
					p = Popen(["glish",'-l',malt.sd+'ldata_malt90_ql.g','-plain',filename])
					p.wait()
				if not quicklook:
					redlog.set_val(ftemp.replace('.rpf',''),"ldata",malt.vnum)
		

def do_gridzilla(source,filenames,lines,freqs,ifs,force=False,quicklook=False):
	# print(over_gzilla_both)
	make_dirs("gridzilla",lines)
	redlog = ReduceLog.ReduceLog()
	if len(filenames) == 2:
		fn1 = filenames[0]
		fn2 = filenames[1]
		gzilla_do1 = redlog.check_val(fn1.replace(".sdfits",""),"gzilla",malt.vnum) 
		gzilla_do2 = redlog.check_val(fn2.replace(".sdfits",""),"gzilla",malt.vnum) 
	for fn in filenames:
		fail_flag = False
		gzilla_needed = redlog.check_val(fn.replace(".sdfits",""),"gzilla",malt.vnum) 
		for i,line,freq in zip(ifs,lines,freqs):
			filein  = fn.replace(".sdfits","")+'_'+line+'.sdfits'
			fileout = fn.replace(".sdfits","")+'_'+line
			#print(filein)
			if gzilla_needed or force:
				q = Popen(["glish",'-l',malt.sd+'gzill_malt90.g',\
					line,str(freq),str(i-1),fileout,filein])
				q.wait() 
				try:
					shutil.move(malt.data_dir+'/gridzilla/'+line+'/'+fileout+'.beamRSS.fits',malt.data_dir+'/gridzilla/'+line+'/'+'flotsam')
					shutil.move(malt.data_dir+'/gridzilla/'+line+'/'+fileout+'.spectracounts.fits',malt.data_dir+'/gridzilla/'+line+'/'+'flotsam')
				except IOError:
					fail_flag = True
					print("Failed to move flotsam files")
		if (not fail_flag) and (not quicklook):
			redlog.set_val(fn.replace(".sdfits",""),"gzilla",malt.vnum)
	if len(filenames) == 2:
		for i,line,freq in zip(ifs,lines,freqs):
			file1   = fn1.replace(".sdfits","")+'_'+line+'.sdfits'
			file2   = fn2.replace(".sdfits","")+'_'+line+'.sdfits'
			fileout = source+'_'+line
			if gzilla_do1 or gzilla_do2 or force:
				q = Popen(["glish",'-l',malt.sd+'gzill_malt90.g',\
					line,str(freq),str(i-1),fileout,file1,file2])
				q.wait()
				try:
					shutil.move(malt.data_dir+'/gridzilla/'+line+'/'+fileout+'.beamRSS.fits',malt.data_dir+'/gridzilla/'+line+'/'+'flotsam')
					shutil.move(malt.data_dir+'/gridzilla/'+line+'/'+fileout+'.spectracounts.fits',malt.data_dir+'/gridzilla/'+line+'/'+'flotsam')
				except IOError:
					print("Failed to move flotsam files")
					#This does not return a useful error
def create_source_folder(source,lines,force=False,quicklook=False):
	"""Create a folder of symlinks for each source"""
	redlog = ReduceLog.ReduceLog()
	try:
		os.mkdir(malt.data_dir+"sources")
	except OSError:
		pass
	for line in lines:
		try:
			os.mkdir(malt.data_dir+'sources/'+source)
		except OSError:
			pass
		filename = source+'_'+line+'_MEAN.fits'
		src  = malt.data_dir+'gridzilla/'+line+'/'+filename
		targ = malt.data_dir+'sources/'+source+'/'+filename
		try:
			os.symlink(src,targ)
		except OSError:
			pass
		hdulist = pyfits.open(targ,mode='update')
		prihdr = hdulist[0].header
		prihdr.update('M90PIPEV',malt.vnum,'Malt90 Pipeline Version')
		prihdr.update('BUNIT','K','Antenna Temperature')
		hdulist.flush()
		#except:
			#print("Header update failed")
	files_involved = [source+"_GLat",source+"_GLon"]
	for file_involved in files_involved:
		#print(file_involved)
		try:
			if not quicklook:
				redlog.set_val(file_involved,"arrange",malt.vnum)
		except:
			print("Failed to update log for "+file_involved+". Maybe it does not exist?")
	
def do_mommaps(source,filenames,lines,force=False,quicklook=False):
	"""Make moment maps for source and place them correctly
	Determine a source velocity (currently from table)
	Make moment map for GLon, GLat, combined in mommaps/line
	Copy the combined moment map into sources/ folder
	"""

	redlog = ReduceLog.ReduceLog()

	for file_involved in filenames:
	       	mommap_needed = False
		direction = file_involved.partition('_')[2]
		print(direction)
		if mommap_needed == False:
			mommap_needed = redlog.check_val(file_involved,"mommaps",malt.vnum) 
		if mommap_needed or force:
			print("I am doing a moment map")
			if quicklook:
				moment_map.do_source(source,lines,direction=direction,auto=True)
			else:
				moment_map.do_source(source,lines,direction=direction)
		for line in lines:
			momsrc = malt.data_dir+'mommaps/'+line+'/'+source+'_'+line+'_mommaps'
			momtarg = malt.data_dir+'sources/'+source+'/'+source+'_'+line+'_mommaps'
			try:
				shutil.copytree(momsrc,momtarg)
			except OSError:
				print("Failed to copy moment maps")
				pass
		if not quicklook:
			redlog.set_val(file_involved,"mommaps",malt.vnum)

if __name__ == '__main__':
	main()

