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
	       Valid options: ldata,gzilla,mommaps,reorg
-j : JustIm -- just remake the integrated intensity map verification image, not
               the spectra. A single-use option to remake year1 and year2 images
	       after the relevant files have already been moved out.
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
import preprocess_malt
import datetime

def main():
	try:
		opts,args=getopt.getopt(sys.argv[1:], "n:s:af:i:jh")
	except getopt.GetoptError,err:
		print str(err)
		print(__doc__)
		sys.exit(2)
	force_list  = []
	ignore_list = []
	do_source = False
	do_date   = False
	do_all    = False
	justim    = False

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
			print("Reducing a night of data -- "+date)
		elif o == "-s":
			do_source = True
			source = a
			print("Reducing a single source")
		elif o == "-a":
			#I think force does not work well with do_all
			do_all = True
			print("Reducing all undone sources")
		elif o == "-j":
			justim = True
		elif o == "-h":
			print(__doc__)
			sys.exit(1)
		else:
			assert False, "unhandled option"
			print(__doc__)
			sys.exit(2)

	#Make sure that log is up-to-date.
	#This step is quick, so just do +/- on current UT. 
	now_utc = datetime.datetime.utcnow()
	today = datetime.date(now_utc.year,now_utc.month,now_utc.day)
	plus1 = today + datetime.timedelta(days=1)
	minu1 = today - datetime.timedelta(days=1)

	for dates in [today,plus1,minu1]:
		files_to_process = preprocess_malt.get_new_files(
			           dates.isoformat(),in_middle_of_obs = True)
		preprocess_malt.rename_files(files_to_process)

	redlog = ReduceLog.ReduceLog()
	if do_date:
		sources = redlog.find_files_with_date(date)
	elif do_source:
		sources = [source]
	elif do_all:
		if force_list:
			force = force_list[0]
		else:
			force = None
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
			do_reduction(one_source,force_list,ignore_list,justim=justim)
		except (IOError,OSError),e:
			print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
			print("OS Error while processing "+one_source)
			print(e)


def do_reduction(source,force_list=[],ignore_list=[],
		 quicklook=False,onlyone=None,justim=False):
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
	### Always do this step unless explicitly ignored###
	if 'reorg' in ignore_list:
		pass
	else:
		create_source_folder(source,lines,quicklook=quicklook)

	### Do moment maps ###
	if onlyone:
		filenames = [source+'_'+onlyone]
	else:
		filenames = [source+"_GLat",source+"_GLon",source]
	if 'mommaps' in force_list:
		#pass #Temporarily disable. Restore this!
                try:
			do_mommaps(source,filenames,lines,force=True,
			   quicklook=quicklook)
		except IOError:
			pass
	elif 'mommaps' in ignore_list:
		pass
	else:
		try:
			do_mommaps(source,filenames,lines,force=False,
				   quicklook=quicklook)
		except IOError:
			pass
	### Make Verification Images ###
	### Always do this step too  ###
	make_verification_plots(source,direction=onlyone,justim=justim)

def setup_lines(quicklook=False,patricio_flag=False):
	if quicklook:
		### Lines to use for quick look ###
		lines = ["hnc","hcop"]
		freqs = [90663.572,89188.526]
		ifs   = [7,9]
	elif patricio_flag:
        	### Special code for Patricio ###
		lines = ["n2hp","13cs","h41a","ch3cn",
			 "hc3n","13c34s","hnc","hc13ccn",
			 "hcop","hcn","hnco413","hnco404",
			 "c2h","hn13c"]
		freqs = [96741.420, 96412.961, 95914.310, 95169.516, 
			 94410.895, 93870.098, 93173.772, 92494.303, 
			 91985.316, 90978.989, 90926.036, 90663.572, 
			 90593.059, 89188.526]	
		ifs   = [1,2,3,4,5,6,7,8,9,10,11,12,13,14]
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
	print("Entering Livedata function")
	make_dirs("livedata",lines)
	redlog = ReduceLog.ReduceLog()
	for filename in filenames:
		print(filename)
		print(malt.data_dir+'renamed/'+filename)
		if os.path.exists(malt.data_dir+'renamed/'+filename):
			print("File Exists")
			ftemp = filename
			ldata_needed = redlog.check_val(ftemp.replace(
					'.rpf',''),"ldata",malt.vnum["ldata"])
			print(force)
			if ldata_needed or force:
				#Could use check_call here to see if this dies
				if not quicklook:
					print("Doing Livedata")
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
				except OSError:
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
				except OSError:
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
		#src  = '../../gridzilla/'+line+'/'+filename
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
		direction = file_involved.partition('_')[2].lstrip('_2').lstrip('_3')
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
						     direction=direction,
						     auto=False)
		for line in lines:
			if quicklook:
				endpart = source+'_'+direction+'_'+line+'_mommaps'
			else:
				endpart = source+'_'+line+'_mommaps'
			
			momsrc = malt.data_dir+'mommaps/'+line+'/'+endpart
			momtarg = malt.data_dir+'sources/'+source+'/'+endpart
			try:
				shutil.rmtree(momtarg)
			except OSError:
				print("No moment map directory to remove.")
				print("Did not remove: "+momtarg)
			try:
				shutil.copytree(momsrc,momtarg)
			except OSError:
				print("Failed to copy moment maps")
				print("From "+momsrc)
				print("To "+momtarg)
		if not quicklook:
			redlog.set_val(file_involved,"mommaps",
				       malt.vnum["mommaps"])
def robust_getdata(base):
	fail = False
	data = 0
	h = 0
	print("Entering robust_getdata")
	try:
		data,h = pyfits.getdata(malt.data_dir+base,header=True)
		#print("Trying base dir...")
	except:
		try:
			#print("Trying year1 dir...")
			data,h = pyfits.getdata(malt.data_dir_y1+base,header=True)
		except:
			try:
				data,h = pyfits.getdata(malt.data_dir_y2+base,header=True)
			except:
				fail = True
	return(data,h,fail)

def make_verification_plots(source,direction=None,justim=False):
	"""Make an image of the 0th moment map
	and of the sepctrum at the peak position with
	a line indicating the central velocity"""
	lines = ["hcop","hnc"]

#	print("This is justim")
	#print(justim)
	if not direction:
		direction = ""
	else:
		direction = "_"+direction
	#print("This is direction")
	#print(direction)
	failflag = False
	for line in lines:
		print(source)
		print(line)
		cube,h,fail = robust_getdata('gridzilla/'+line+'/'+source+direction+"_"+line+"_MEAN.fits")
				
		if fail:
			failflag = True
			break
		snr = np.nan_to_num(robust_getdata('mommaps/'
						   +line+'/'+source+direction
						   +"_"+line+"_mommaps/"
						   +source+direction+"_"
						   +line+"_snr0.fits")[0])
		mask = np.zeros(snr.shape)
		try:
			mask[3:28,3:28] = 1
		except IndexError:
			failflag = True
			break
		d,hmom,fail = robust_getdata('mommaps/'+line
					+'/'+source+direction+"_"+line
					+"_mommaps/"+source+direction+"_"
					+line+"_mom0.fits")
		d = np.nan_to_num(d)
		plt.clf()
		im = idl_stats.blur_image(d*mask,3)
		plt.imshow(im[::-1,:])#, vmin=0, vmax = 19)
#		print(im[::-1,:])
		a = plt.colorbar()
		a.set_label("K km/s")
		plt.title(source+" "+direction+" "+line
			  +" integrated intensity")
		plt.xticks((0,5,10,15,20,25),(-13,-8,-3,2,7,12))
		plt.yticks((0,5,10,15,20,25),(-13,-8,-3,2,7,12))
		plt.ylabel("Galactic Latitude Offset [9 arcsec pixels]")
		plt.xlabel("Galactic Longitude Offset [9 arcsec pixels]")
		plt.savefig(malt.data_dir+'verification/'+source+"_"
			    +line+"_mom0"+".png")

		if not justim:
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
				    +line+"_velcheck"+".png")
	if failflag:
		print("=================================================")
		print("Something failed in the reduction of this source.")
		print("Source name: "+source)
		print("The most likely cause is that the data was not   ")
		print("yet synced over from the telescope. Data syncs at")
		print("3, 18, 33, and 48 minutes past the hour. Wait    ")
		print("until a minute past one of these marks and try   ")
		print("your reduction command again.                    ")
		print("=================================================")
	else:
		print("=================================================")
		print("This source ("+source+") appears to have         ")
		print("reduced successfully.                            ")
		print("=================================================")
if __name__ == '__main__':
	main()

