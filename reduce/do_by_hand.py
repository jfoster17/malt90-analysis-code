#!/usr/bin/env python
# encoding: utf-8
"""
do_by_hand.py does special reductions of sources outside
of the pipeline. Typical usage is to combine multiple 
different maps of one source. Files are stored in byhand/
directory to avoid any conflicts in names while keeeping
a simple naming system. 

Options
-h : Display this help
-s : Specify the source name
-v : Specify a central velocity (normally use this)
-f : List raw data files to combine
-p : Do reduction for Patricio
"""

import sys,os,getopt
import pyfits
import numpy as np
import shutil
import subprocess as sp
import malt_params as malt
import idl_stats

import preprocess_malt,reduce_malt,ReduceLog,moment_map

def main():
	try:
		opts,args = getopt.getopt(sys.argv[1:],"s:f:v:hp")
	except getopt.GetoptError,err:
		print(str(err))
		print(__doc__)
		sys.exit(2)
	patricio_flag = False
	velocity = -999
	for o,a in opts:
		if o == "-p":
			patricio_flag = True
		if o == "-s":
			source = a
		if o == "-f":
			input_files = a.split(",")
			print(input_files)
		if o == "-h":
			print(__doc__)
			sys.exit(1)
		if o == "-v":
			velocity = a
	working_names = copy_files(source,input_files)
	reduce_map(source,working_names,patricio_flag,velocity)

def copy_files(source,input_files):
#	print(source)
	working_names = []
	for i,input_file in enumerate(input_files):
		working_name = source+"_"+str(i)
		working_names.append(working_name)
	#	print(malt.source_dir+input_file)
		#print(malt.byhand_rename_dir+working_name+".rpf")
		shutil.copyfile(malt.source_dir+input_file,
				malt.byhand_rename_dir+working_name+".rpf")
	return(working_names)


def do_gridzilla(source,working_names,patricio_flag=False):
	### Do Gridzilla ###
	lines,freqs,ifs = reduce_malt.setup_lines(patricio_flag=patricio_flag)

	# First do all the individual files so we 
	# can check them later for quality
	for fname in working_names:
		for i,line,freq in zip(ifs,lines,freqs):
			fout = fname+"_"+line
			fin  = fname+"_"+line+".sdfits"
			q = sp.Popen(["glish",'-l',malt.sd+
				   'gzill_malt90_byhand.g',line,
				   str(freq),str(i-1),fout,
				   fin])
			q.wait()

	#Then do all the files together
	#Terrible hack to make a flexible
	#number of input files work with glish
	for i,line,freq in zip(ifs,lines,freqs):
		fileout = source+'_'+line
		af = []
		for fname in working_names:
			af.append(fname+"_"+line+".sdfits")

		print(af)
		nfiles = len(af)
		if nfiles == 1:
			q = sp.Popen(["glish",'-l',malt.sd+
				      'gzill_malt90_byhand.g',
				      line, str(freq),str(i-1),
				      fileout,af[0]])
		elif nfiles == 2:
			q = sp.Popen(["glish",'-l',malt.sd+
				      'gzill_malt90_byhand.g',
				      line, str(freq),str(i-1),
				      fileout,af[0],af[1]])
		elif nfiles == 3:
			q = sp.Popen(["glish",'-l',malt.sd+
				      'gzill_malt90_byhand.g',
				      line, str(freq),str(i-1),
				      fileout,af[0],af[1],af[2]])
		elif nfiles == 4:
			q = sp.Popen(["glish",'-l',malt.sd+
				      'gzill_malt90_byhand.g',
				      line, str(freq),str(i-1),
				      fileout,af[0],af[1],af[2],af[3]])
		elif nfiles == 5:
			q = sp.Popen(["glish",'-l',malt.sd+
				      'gzill_malt90_byhand.g',
				      line, str(freq),str(i-1),
				      fileout,af[0],af[1],af[2],
				      af[3],af[4]])
		elif nfiles == 6:
			q = sp.Popen(["glish",'-l',malt.sd+
				      'gzill_malt90_byhand.g',
				      line, str(freq),str(i-1),
				      fileout,af[0],af[1],af[2],
				      af[3],af[4],af[5]])

		q.wait()


def reduce_map(source,working_names,patricio_flag=False,velocity=-999):

	print(velocity)
	### Do Livedata ###
	for filename in working_names:
		p = sp.Popen(["glish",'-l', malt.sd+'ldata_malt90_byhand.g',
                          '-plain',filename+".rpf"])
		p.wait()

	do_gridzilla(source,working_names,patricio_flag = patricio_flag)

	lines,freqs,ifs = reduce_malt.setup_lines(patricio_flag=patricio_flag)
	### Do Source Folder ###
	for line in lines:
		try:
			os.mkdir(malt.byhand_data_dir+
				 'sources/'+source)
		except OSError:
			pass
		filename = source+'_'+line+'_MEAN.fits'
		src  = malt.byhand_data_dir+'gridzilla/'+line+'/'+filename
		targ = malt.byhand_data_dir+'sources/'+source+'/'+filename
		try:
			os.unlink(targ)
		except OSError:
			pass
		try:
			os.symlink(src,targ)
		except OSError:
			pass
		hdulist = pyfits.open(targ,mode='update')
		prihdr = hdulist[0].header
		prihdr.update('M90PIPEV',malt.vnum["gzilla"],
			      'Malt90 Pipeline Version')
		prihdr.update('BUNIT','K',
			      'Antenna Temperature')
		prihdr.update('BMAJ',0.0105,
			      'Beam major FWHM (degrees)')
		prihdr.update('BMIN',0.0105,
			      'Beam minor FWHM (degrees)')
		hdulist.flush()
	if velocity != -999:
		moment_map.do_source(source,lines,
			     auto=False,altdir=malt.byhand_data_dir,vel=velocity)
	else:
		moment_map.do_source(source,lines,
			     auto=True,altdir=malt.byhand_data_dir)
	for line in lines:
		endpart = source+'_'+line+'_mommaps'
		
		momsrc = malt.byhand_data_dir+'mommaps/'+line+'/'+endpart
		momtarg = malt.byhand_data_dir+'sources/'+source+'/'+endpart
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


if __name__ == '__main__':
	main()

