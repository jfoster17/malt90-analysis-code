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

base_data = "/DATA/MALT_1/MALT90/data/"

def get_velocity(source):
	"""Get a velocity for a source.
	For now, use tabulated value.
	Later, this function will find a velocity.
	"""
	velocity = 30.6 #Dummy for G303.930
	return(velocity)
	
def do_source(source,lines):
	print("Sourcename:")
	print(source)
	central_velocity = get_velocity(source)
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
	
	
def get_filename(source,line):
	filename = source+"_"+line+"_MEAN.fits"
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
	if central_velocity: #This signals a good velocity
		hispec = np.nonzero(vel >= central_velocity)[0]
		cenchan = hispec[0]
		minchan = max([cenchan - n_pad,n_edge])
		maxchan = min([cenchan + n_pad,nchan-n_edge])
	else:
		pass
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

def main():
	pass

if __name__ == '__main__':
	main()

