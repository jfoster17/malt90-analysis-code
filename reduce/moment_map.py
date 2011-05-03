#!/usr/bin/env python
# encoding: utf-8
"""
moment_map.py
"""

import sys
import os
import pyfits
import numpy as np
import idl_stats
import pylab
import numpy.ma as ma
import reduce_malt
import scipy.ndimage

base_data = "/DATA/MALT_1/MALT90/data/"

def get_velocity(source,auto=True,direction=None):
	"""Get a velocity for a source.
	For now, use tabulated value.
	Later, this function will find a velocity.
	"""
	if auto:
		velocity = identify_velocity(source,direction=direction)
	else:
		path_to_vel = os.path.join(reduce_malt.sd,'malt90_velocities_year1.txt')
		f = open(path_to_vel,'r')
		for line in f:
			if line.split()[0].strip() == source:
				velocity = float(line.split()[1])

		f.close()
		#print(velocity)
	return(velocity)

def identify_velocity(source,minchan = 200,maxchan = 3896,sig=5,direction=None):
	"""Identify a source velocity based on HCO+
	Later I want to add HNC
	"""
	try:
		infile = get_filename(source,"hcop",direction=direction)
		d,h = pyfits.getdata(infile,header=True)
	except OSError:
		print("Failed to open datacube "+infile)
		return(0)
	nglat = d.shape[1]
	nglon = d.shape[2]
	nspec = d.shape[0]
	vel = (np.arange(nspec)+1-h['CRPIX3'])*h['CDELT3']+h['CRVAL3']
	sf = 11
	threshold = 4
	channel_trim = 7
	edge = 3
	max_chan = np.zeros((nglat,nglon))
	for x in range(edge,nglat-edge):
		for y in range(edge,nglon-edge):
			spec = d[minchan:maxchan,x,y]
			smoothspec = smooth(spec,window_len=sf,window='hamming')
			mean,sigma = idl_stats.iterstat(smoothspec)
			goodsignal = np.where(smoothspec > threshold*sigma,1,0)
			goodsignal = scipy.ndimage.binary_erosion(goodsignal,structure=np.ones(channel_trim))
			maskedsignal = goodsignal*smoothspec
			max_chan[x,y] = np.argmax(maskedsignal)
	max_chan = np.extract(max_chan > 0,max_chan)
	best_chan = int(np.median(max_chan))+minchan
	return(vel[best_chang]/1000.)

def do_source(source,lines,direction=None,auto=False):
	print("Sourcename:")
	print(source)
	central_velocity = get_velocity(source,auto,direction)
	create_basic_directories(lines)
	#create_output_directories(source,lines)
	for line in lines:
		infile = get_filename(source,line)
		print(infile)
		out_base = infile[:-9].replace("gridzilla","mommaps")
		out_dir = source+"_"+line+"_mommaps"
		try:
			output_dir = os.path.join(base_data,"mommaps",line,out_dir)
			os.mkdir(output_dir)
		except OSError:
			pass
		make_moment_maps(infile,out_base,output_dir,central_velocity=central_velocity*1000)

def create_basic_directories(lines):
	"""Create subdirectories under mommaps"""
	for line in lines:
		try:
			moment_dir = os.path.join(base_data,"mommaps",line)
			os.mkdir(moment_dir)
		except OSError:
			pass
	
	
def get_filename(source,line,direction=None):
	if not direction:
		filename = source+"_"+line+"_MEAN.fits"
	else:
		filename = source+"_"+direction+"_"+line+"_MEAN.fits"
	full_path = os.path.join(base_data,"gridzilla",line,filename)
	return(full_path)

def calculate_moments(d,minchan=False,maxchan=False,vel=False,bestmask=False,mask=False):
	"""This function actually calculates moments"""
	nglat = d.shape[1]
	nglon = d.shape[2]
	nspec = d.shape[0]
	maps = np.zeros((nglat,nglon),dtype={'names':['mean','sd','errmn',
		'errsd','skew','kurt','error','intint','npix'],
		'formats':['f4','f4','f4','f4','f4','f4','f4','f4','f4']})
	#These definitions for mask seem backward but are correct.
	noise_portion = ma.masked_where(mask == 1,d)
	good_d = d[minchan:maxchan,...]
	mask2 = mask[minchan:maxchan,...]
	print(minchan)
	print(maxchan)
	signal_portion = ma.masked_where(mask2 == 0,good_d)
	maps['error']  = ma.std(noise_portion,axis=0)
	maps['intint'] = ma.sum(signal_portion,axis=0)
	for x in range(nglat):
		for y in range(nglon):
			fullspec = d[...,x,y]#Exract a single spectrum
			ind = np.arange(nspec)
			velmask = mask[minchan:maxchan,x,y]
			if np.sum(velmask) != 0:
				velmask = bestmask
				npix = max(np.sum(velmask),1)
			ind = ind[velmask > 0]
			sigma = maps['error'][x,y]
			if ind.size > 2 and (sigma > 0):
				mom = idl_stats.wt_moment(vel[ind],fullspec[ind],
						errors = np.zeros(ind.size)+sigma)
				maps['mean'][x,y]  = mom['mean']
				maps['sd'][x,y]    = mom['stdev']
				maps['errmn'][x,y] = mom['errmn']
				maps['errsd'][x,y] = mom['errsd']
				maps['npix'][x,y]  = npix
			else:
				maps['mean'][x,y]  = np.nan
				maps['sd'][x,y]    = np.nan
				maps['errmn'][x,y] = np.nan
				maps['errsd'][x,y] = np.nan
				maps['npix'][x,y]  = np.nan
	return(maps)

def make_moment_maps(infile,out_base,output_dir,central_velocity=False,second=False):
	"""Wrapper function to deal with headers
	Call function to make maps
	Save maps
	Central velocity is in m/s
	n_pad ~ 45 km/s, appropriate for N2H+
	"""
	n_edge = 10 #Bad/noisy edge chanels
	print("Processing..."+infile)
	
	d,h = pyfits.getdata(infile,header=True)
	nchan = h['NAXIS3']
	vel = (np.arange(nchan)+1-h['CRPIX3'])*h['CDELT3']+h['CRVAL3']
	
	hdout = h
	vwidth = h['CDELT3']
	
	del hdout['CRVAL3']
	del hdout['CRPIX3']
	del hdout['CDELT3']
	del hdout['CTYPE3']
	del hdout['NAXIS3']

	minchan = n_edge
	maxchan = nchan - n_edge 
	
	print("Doing Pre-determine velocity integration...")
	#maps = do_predetermined_velocity(central_velocity,vel,hdout,n_edge,nchan,d,n_pad = 125)
	#save_maps(maps,hdout,out_base,out_dir,vel,minchan,maxchan,vwidth,"fullvel")
	
	maps = do_predetermined_velocity(central_velocity,vel,hdout,n_edge,nchan,d,n_pad = 75)
	print("Maps Made")
	print(out_base)
	print(output_dir)
#	print(hdout)
	save_maps(maps,hdout,out_base,output_dir,vel,minchan,maxchan,vwidth,"medvel")

	#maps = do_predetermined_velocity(central_velocity,vel,hdout,n_edge,nchan,d,n_pad = 25)
	#save_maps(maps,hdout,out_base,out_dir,vel,minchan,maxchan,vwidth,"smallvel")

def do_predetermined_velocity(central_velocity,vel,hdout,n_edge,nchan,d,n_pad):
	"""Make a moment map once we know central velocity."""

	hispec = np.nonzero(vel >= central_velocity)[0]
	cenchan = hispec[0]
	minchan = max([cenchan - n_pad,n_edge])
	maxchan = min([cenchan + n_pad,nchan-n_edge])

	vmin = vel[minchan]/1e3
	vmax = vel[maxchan]/1e3
	print("Velocity Integration Limit: "+str(vmin)+' to '+str(vmax))
	hdout.update('VMIN',vmin,'KM/S')
	hdout.update('VMAX',vmax,'KM/S')
	mask = np.zeros(d.shape,dtype=np.int)
	mask[minchan:maxchan,...] = 1
	bestmask = mask[...,0,0]
	maps = calculate_moments(d,minchan,maxchan,vel,bestmask=bestmask,mask=mask)
	return(maps)

def save_maps(maps,hdout,out_base,out_dir,vel,minchan,maxchan,vwidth,name_mod):
	
	(head,tail) = os.path.split(out_base)
	out_base2 = os.path.join(out_dir,tail)
#	out_temp = os.path.join(out_base2,tail)
	print("Base for output")
	print(out_base2)
#	print(out_temp)
	out_base = out_base2
	badind = np.where((maps['errmn'] > 1e6) | (maps['errsd'] > 1e6)) #This trims out sources with sigma_v > 1000 km/s
	try:
		maps['mean'][badind] = np.nan
		maps['sd'][badind]   = np.nan
	except IndexError:
		pass
	maps['skew'] = maps['skew']/maps['sd']**3
	maps['kurt'] = maps['kurt']/maps['sd']**4 - 3 #Really?
	hdout['NAXIS'] = 2
	hdout.update('CDELT1',hdout['CDELT1'],'DEGREES')
	hdout.update('CDELT2',hdout['CDELT2'],'DEGREES')
	hdout.update('CRVAL1',hdout['CRVAL1'],'DEGREES')
	hdout.update('CRVAL2',hdout['CRVAL2'],'DEGREES')
	hdout.update('BUNIT','KM/S')
	
	maps['intint'] = maps['intint']*vwidth/1e3
	pyfits.writeto(out_base+'mom1'+'.fits',maps['mean']/1e3,hdout,clobber=True)
	pyfits.writeto(out_base+'mom2'+'.fits',maps['sd']/1e3,hdout,clobber=True)
	pyfits.writeto(out_base+'err1'+'.fits',maps['errmn']/1e3,hdout,clobber=True)
	pyfits.writeto(out_base+'err2'+'.fits',maps['errsd']/1e3,hdout,clobber=True)
	hdout.update('BUNIT','NONE')

	pyfits.writeto(out_base+'npix'+'.fits',maps['npix'],clobber=True)
	pyfits.writeto(out_base+'snr0'+'.fits',maps['intint']/(maps['error']*
		np.sqrt(maps['npix'])*vwidth/1e3),hdout,clobber=True)
	
	hdout.update('BUNIT','K.KM/S')
	
	pyfits.writeto(out_base+'mom0'+'.fits',maps['intint'],hdout,clobber=True)
	pyfits.writeto(out_base+'err0'+'.fits',maps['error']*np.sqrt(maps['npix'])*
					vwidth/1e3,hdout,clobber=True)
	pyfits.writeto(out_base+'emap'+'.fits',maps['error']*
					vwidth/1e3,hdout,clobber=True)

def smooth(x,window_len=11,window='hanning'):
	"""smooth the data using a window with requested size.

	This method is based on the convolution of a scaled window with the signal.
	The signal is prepared by introducing reflected copies of the signal 
	(with the window size) in both ends so that transient parts are minimized
	in the begining and end part of the output signal."""

	if x.ndim != 1:
		raise ValueError, "smooth only accepts 1 dimension arrays."

	if x.size < window_len:
		raise ValueError, "Input vector needs to be bigger than window size."
	
	if window_len<3:
		return x

	if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
		raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
	
	
	s=np.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
						#print(len(s))
	if window == 'flat': #moving average
		w=ones(window_len,'d')
	else:
		w=eval('np.'+window+'(window_len)')
		
	y=np.convolve(w/w.sum(),s,mode='same')
	return y[window_len-1:-window_len+1]
	


def main():
	pass

if __name__ == '__main__':
	main()

