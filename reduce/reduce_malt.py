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
-n : Night  -- reduce all files taken on a give date (YYYY-MM-DD)
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
-h : Help   -- display this help

--- Changelog ---
Version 1.6 uses minimal gaussian smoothing in Gridzilla.
Version 1.5 correctly identifies SiO and HN13C lines data 
(they were swapped in previous versions)
Version 1.4 removes the ASAP smoothing (which distorted the 
velocity axis) and uses Livedata to smooth
"""

import sys,os,shutil,getopt
from subprocess import *
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import pyfits
import numpy as np
import ReduceLog
import moment_map
import idl_stats
import malt_params as malt

def main():
	try:
		opts,args=getopt.getopt(sys.argv[1:], "n:s:af:i:h")
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
			print("Forcing reduction of: "+str(force_list))
		elif o == "-i":
			ignore_list = a.split(',')
			print("Not doing: "+str(ignore_list))
		elif o == "-n":
			do_date = True
			date = a
			print("Reducing a night of data -- "+data)
		elif o == "-s":
			do_source = True
			source = a
			print("Reducing a single source")
		elif o == "-a":
			#I think force does not work well with do_all
			do_all = True
			print("Reducing all undone sources")
		elif o == "-h":
			print(__doc__)
			sys.exit(1)
		else:
			assert False, "unhandled option"
			print(__doc__)
			sys.exit(2)

	redlog = ReduceLog.ReduceLog()
	if do_date:
		sources = redlog.find_files_with_date(date)
	elif do_source:
		sources = [source]
	elif do_all:
		if force_list:
			force = force_list[0]
		sources = redlog.find_undone(malt.vnum,force=force)
		sources = list(set(sources)) #Eliminate dupes
		print("Reducing the following source(s): "+str(sources))
	try:
		print("Reducing the following source(s): "+str(sources))
	except UnboundLocalError:
		print("@@@ No sources to reduce @@@")
		print(__doc__)
		sys.exit(2)
	for one_source in sources:
		try:
			do_reduction(one_source,force_list,ignore_list)
		except IOError,e:
			print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
			print("IO Error while processing "+one_source)
			print(e)


def do_reduction(source,force_list=[],ignore_list=[],
		 quicklook=False,onlyone=None):
	"""Main script calling Livedate,Gridzilla,Organize,MomentMaps"""

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
		do_gridzilla(source,filenames,lines,freqs,ifs,force=True,
			     quicklook=quicklook)
	elif 'gzilla' in ignore_list:
		pass
	else:
		do_gridzilla(source,filenames,lines,freqs,ifs,force=False,
			     quicklook=quicklook)
	### Do Reorganization ###
	### Always do this step ###
	create_source_folder(source,lines,quicklook=quicklook)

	### Do moment maps ###
	if onlyone:
		filenames = [source+'_'+onlyone]
	else:
		filenames = [source+"_GLat",source+"_GLon",source]
	if 'mommaps' in force_list:
		do_mommaps(source,filenames,lines,force=True,
			   quicklook=quicklook)
	elif 'mommaps' in ignore_list:
		pass
	else:
		do_mommaps(source,filenames,lines,force=False,
			   quicklook=quicklook)
	### Make Verification Images ###
	### Always do this step too  ###
	make_verification_plots(source)

def setup_lines(quicklook=False):
	if quicklook:
		### Lines to use for quick look ###
		lines = ["hnc","hcop"]
		freqs = [90663.572,89188.526]
		ifs   = [7,9]
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
	### Malt90 backward sources (first day) ###
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
		if dirname=="gridzilla":
			try:
				os.mkdir(malt.data_dir+dirname+"/"+line+"/flotsam")
			except OSError:
				pass
def do_livedata(filenames,lines,force=False,quicklook=False):
	"""Do Livedata if we need to"""
	make_dirs("livedata",lines)
	redlog = ReduceLog.ReduceLog()
	for filename in filenames:
		print(filename)
		if os.path.exists(malt.data_dir+'renamed/'+filename):
			ftemp = filename
			ldata_needed = redlog.check_val(ftemp.replace(
					'.rpf',''),"ldata",malt.vnum["ldata"])
			if ldata_needed or force:
				#Could use check_call here to see if this dies
				if not quicklook:
					p = Popen(["glish",'-l',
						   malt.sd+'ldata_malt90.g',
						   '-plain',filename])
					p.wait()
				else:
					p = Popen(["glish",'-l',
						   malt.sd+'ldata_malt90_ql.g',
						   '-plain',filename])
					p.wait()
				if not quicklook:
					redlog.set_val(ftemp.replace('.rpf',
						       ''),"ldata",
						       malt.vnum["ldata"])
		

def do_gridzilla(source,filenames,lines,freqs,ifs,force=False,quicklook=False):
	"""Do Gridzilla on GLat, GLon, and Both"""
	make_dirs("gridzilla",lines)
	redlog = ReduceLog.ReduceLog()
	if len(filenames) == 2:
		fn1 = filenames[0]
		fn2 = filenames[1]
		gzilla_do1 = redlog.check_val(fn1.replace(".sdfits",""),
					      "gzilla",malt.vnum["gzilla"]) 
		gzilla_do2 = redlog.check_val(fn2.replace(".sdfits",""),
					      "gzilla",malt.vnum["gzilla"]) 
	for fn in filenames:
		fail_flag = False
		gzilla_needed = redlog.check_val(fn.replace(".sdfits",""),
						 "gzilla",malt.vnum["gzilla"]) 
		for i,line,freq in zip(ifs,lines,freqs):
			filein  = fn.replace(".sdfits","")+'_'+line+'.sdfits'
			fileout = fn.replace(".sdfits","")+'_'+line
			if gzilla_needed or force:
				q = Popen(["glish",'-l',
					   malt.sd+'gzill_malt90.g',line,
					   str(freq),str(i-1),fileout,filein])
				q.wait() 
				try:
					beam = fileout+'.beamRSS.fits'
					counts = fileout+'.spectracounts.fits'
					os.rename(malt.data_dir
						    +'/gridzilla/'+line+'/'
						    +beam,
						    malt.data_dir+'/gridzilla/'
						    +line+'/flotsam/'+beam)
					os.rename(malt.data_dir
						    +'/gridzilla/'+line+'/'
						    +counts,
						    malt.data_dir
						    +'/gridzilla/'+line
						    +'/flotsam/'+counts)
				except IOError:
					fail_flag = True
					print("Failed to move flotsam files")
		if (not fail_flag) and (not quicklook):
			redlog.set_val(fn.replace(".sdfits",""),
				       "gzilla",malt.vnum["gzilla"])
	if len(filenames) == 2:
		for i,line,freq in zip(ifs,lines,freqs):
			file1   = fn1.replace(".sdfits","")+'_'+line+'.sdfits'
			file2   = fn2.replace(".sdfits","")+'_'+line+'.sdfits'
			fileout = source+'_'+line
			if gzilla_do1 or gzilla_do2 or force:
				q = Popen(["glish",'-l',malt.sd
					   +'gzill_malt90.g',line,str(freq),
					   str(i-1),fileout,file1,file2])
				q.wait()
				try:
					beam = fileout+'.beamRSS.fits'
					counts = fileout+'.spectracounts.fits'
					os.rename(malt.data_dir
						    +'/gridzilla/'+line+'/'
						    +beam,
						    malt.data_dir+'/gridzilla/'
						    +line+'/flotsam/'+beam)
					os.rename(malt.data_dir
						    +'/gridzilla/'+line+'/'
						    +counts,
						    malt.data_dir
						    +'/gridzilla/'+line
						    +'/flotsam/'+counts)
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
		#src  = malt.data_dir+'gridzilla/'+line+'/'+filename
		src  = '../../gridzilla/'+line+'/'+filename
		targ = malt.data_dir+'sources/'+source+'/'+filename
		try:
			os.unlink(targ)
		except OSError:
			pass
		try:
			os.symlink(src,targ)
		except OSError:
			pass
		if not quicklook:
			hdulist = pyfits.open(targ,mode='update')
			prihdr = hdulist[0].header
			prihdr.update('M90PIPEV',malt.vnum["gzilla"],
				      'Malt90 Pipeline Version')
			prihdr.update('BUNIT','K','Antenna Temperature')
			prihdr.update('BMAJ',0.0105,
				      'Beam major FWHM (degrees)')
			prihdr.update('BMIN',0.0105,
				      'Beam minor FWHM (degrees)')
			hdulist.flush()
		#except:
			#print("Header update failed")
	files_involved = [source+"_GLat",source+"_GLon"]
	for file_involved in files_involved:
		#print(file_involved)
		try:
			if not quicklook:
				redlog.set_val(file_involved,"arrange",
					       malt.vnum["arrange"])
		except:
			print("Failed to update log for "+file_involved
			      +". Maybe it does not exist?")
	
def do_mommaps(source,filenames,lines,force=False,quicklook=False):
	"""Make moment maps for source and place them correctly
	Determine a source velocity (currently from table)
	Make moment map for GLon, GLat, combined in mommaps/line
	Copy the combined moment map into sources/ folder
	"""
	print("Starting moment map creation...")
	redlog = ReduceLog.ReduceLog()

	fullmap_needed = False
	for file_involved in filenames:
		part_file = False
		if ('GLat' in file_involved) or ('GLon' in file_involved):
			part_file = True
		mommap_needed = False
		direction = file_involved.partition('_')[2].lstrip('_2')
		print("Direction = "+direction)
		print("Length of direction = ",len(direction)) 
		print("File involved = "+file_involved)
		if mommap_needed == False:
			mommap_needed = redlog.check_val(file_involved,
							 "mommaps",
							 malt.vnum["mommaps"]) 
			
			if part_file:
				fullmap_needed = True
		if mommap_needed or force or (not part_file and fullmap_needed):
			if quicklook:
				moment_map.do_source(source,lines,
						     direction=direction,
						     auto=True)
			else:
				moment_map.do_source(source,lines,
						     direction=direction)
		for line in lines:
			endpart = source+'_'+line+'_mommaps'
			momsrc = malt.data_dir+'mommaps/'+line+'/'+endpart
			momtarg = malt.data_dir+'sources/'+source+'/'+endpart
			try:
				shutil.rmtree(momtarg)
			except OSError:
				print("No moment map directory to remove.")
			try:
				shutil.copytree(momsrc,momtarg)
			except OSError:
				print("Failed to copy moment maps")
				print("From "+momsrc)
				print("To "+momtarg)
		if not quicklook:
			redlog.set_val(file_involved,"mommaps",
				       malt.vnum["mommaps"])

def make_verification_plots(source,direction=None):
	"""Make an image of the 0th moment map
	and of the sepctrum at the peak position with
	a line indicating the central velocity"""
	lines = ["hcop","hnc"]

	if not direction:
		direction = ""
	else:
		direction = "_"+direction
	for line in lines:
		cube,h = pyfits.getdata(malt.data_dir+'gridzilla/'
					+line+'/'+source+direction+"_"
					+line+"_MEAN.fits",header=True)
		snr = np.nan_to_num(pyfits.getdata(malt.data_dir+'mommaps/'
						   +line+'/'+source+direction
						   +"_"+line+"_mommaps/"
						   +source+direction+"_"
						   +line+"_snr0.fits"))
		mask = np.zeros(snr.shape)
		mask[3:28,3:28] = 1
		d,hmom = pyfits.getdata(malt.data_dir+'mommaps/'+line
					+'/'+source+direction+"_"+line
					+"_mommaps/"+source+direction+"_"
					+line+"_mom0.fits",header=True)
		plt.clf()
		plt.imshow(idl_stats.blur_image(d*mask,3))
		a = plt.colorbar()
		a.set_label("K km/s")
		plt.title(source+" "+direction+" "+line
			  +" integrated intensity")
		plt.ylabel("Galactic Latitude Offset [9 arcsec pixels]")
		plt.xlabel("Galactic Longitude Offset [9 arcsec pixels]")
		plt.savefig(malt.data_dir+'verification/'+source+"_"
			    +line+"_mom0"+direction+".png")
	
		nspec = h['NAXIS3']
		vmin = hmom['VMIN']
		vmax = hmom['VMAX']
		vcen = np.average([vmin,vmax])
		vel = ((np.arange(nspec)+1-h['CRPIX3'])*h['CDELT3']
		       +h['CRVAL3'])/1e3
	
		snr_smooth = idl_stats.blur_image(snr,5)
	
		peak_pix = np.argmax(snr_smooth*mask)
		pid = np.unravel_index(peak_pix,snr.shape)

		spectra = idl_stats.smooth(cube[...,pid[0],pid[1]])
	
		plt.clf()
		plt.plot(vel,spectra)
		plt.title(source+" "+direction+" "+line
			  +" at maximum integrated intensity")
		plt.axvline(x=vcen,color='r')
		plt.xlim((vcen-100,vcen+100))
		plt.xlabel("Velocity [km/s]")
		plt.ylabel("T [K]")
		plt.savefig(malt.base+'data/verification/'+source+"_"
			    +line+"_velcheck"+direction+".png")


if __name__ == '__main__':
	main()

