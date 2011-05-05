#!/usr/bin/env python
# encoding: utf-8
"""
quicklook_malt.py is a quick reduction of hco+ and hnc to 
check that we actually see a source.

Call with the name of your source and which direction you want
to look at (leave off direction to combine both)
quicklook_malt.py -s G320.020+00.100 -d GLat

Options
-s : Source The full name of the source you would like to reduce.
-d: Direction The scan direction (GLat or GLon) to reduce. 
    Leave blank to reduce both (reducing both maps takes longer).

Verification images are created in /DATA/MALT_1/MALT90/data/verification/
"""

import sys,os,getopt
import pyfits
import numpy as np

import datetime
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import malt_params as malt
import idl_stats
import preprocess_malt,reduce_malt,ReduceLog

def main():
	try:
		opts,args = getopt.getopt(sys.argv[1:],"s:d:")
	except getopt.GetoptError,err:
		print(str(err))
		print(__doc__)
		sys.exit(2)
	direction = None
	for o,a in opts:
		if o == "-s":
			source = a
		if o == "-d":
			direction = a

	#Make sure that log is up-to-date.
	#This step is quick, so just do +/-
	#on current UT. 
	now_utc = datetime.datetime.utcnow()
	today = datetime.date(now_utc.year,now_utc.month,now_utc.day)
	plus1 = today + datetime.timedelta(days=1)
	minu1 = today - datetime.timedelta(days=1)

	for dates in [today,plus1,minu1]:
		files_to_process = preprocess_malt.get_new_files(dates.isoformat(),in_middle_of_obs = True)
		preprocess_malt.rename_files(files_to_process)
		
	#redlog = ReudceLog.ReduceLog()
	reduce_malt.do_reduction(source,ignore_list = ['ldata','gzilla'],force_list=['mommaps'],quicklook=True,onlyone=direction)
#	reduce_malt.do_reduction(source,force_list=['ldata','gzilla','mommaps'],quicklook=True,onlyone=direction)
	make_verification_plots(source,direction)

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
		cube,h = pyfits.getdata(malt.data_dir+'gridzilla/'+line+'/'+source+direction+"_"+line+"_MEAN.fits",header=True)
		snr = np.nan_to_num(pyfits.getdata(malt.data_dir+'mommaps/'+line+'/'+source+direction+"_"+line+"_mommaps/"+source+direction+"_"+line+"_snr0.fits"))
		mask = np.zeros(snr.shape)
		mask[2:28,2:28] = 1
		d,hmom = pyfits.getdata(malt.data_dir+'mommaps/'+line+'/'+source+direction+"_"+line+"_mommaps/"+source+direction+"_"+line+"_mom0.fits",header=True)
		plt.clf()
		plt.imshow(d*mask)
		a = plt.colorbar()
		a.set_label("K km/s")
		plt.title(source+" "+direction+" "+line+" integrated intenxity")
		plt.ylabel("Galactic Latitude Offset [9 arcsec pixels]")
		plt.xlabel("Galactic Longitude Offset [9 arcsec pixels]")
		plt.savefig(malt.data_dir+'verification/'+source+"_"+line+"_mom0"+direction+".png")
	
		nspec = h['NAXIS3']
		vmin = hmom['VMIN']
		vmax = hmom['VMAX']
		vcen = np.average([vmin,vmax])
		vel = ((np.arange(nspec)+1-h['CRPIX3'])*h['CDELT3']+h['CRVAL3'])/1e3
	
		snr_smooth = idl_stats.blur_image(snr,3)
	
		peak_pix = np.argmax(snr_smooth*mask)
		pid = np.unravel_index(peak_pix,snr.shape)

		#print(pid)
		spectra = idl_stats.smooth(cube[...,pid[0],pid[1]])
		#print(spectra)
		#print(len(vel),len(spectra))
	
		plt.clf()
		plt.plot(vel,spectra)
		plt.title(source+" "+direction+" "+line+" at maximum integrated intensity")
		plt.axvline(x=vcen,color='r')
		plt.xlim((vcen-100,vcen+100))
		plt.xlabel("Velocity [km/s]")
		plt.ylabel("T [K]")
		plt.savefig(malt.base+'data/verification/'+source+"_"+line+"_velcheck"+direction+".png")
		

if __name__ == '__main__':
	main()

