#!/usr/bin/env python
# encoding: utf-8
"""
moment_map.py
"""

import sys,os
import pyfits
import numpy as np
import numpy.ma as ma
import scipy.ndimage
import pylab
import idl_stats
import malt_params as malt

def get_velocity(source,auto=False,direction=None):
	"""Get a velocity for a source
	using either the HCO+ and HNC lines or
	looking up from a table. Try table first
	unless specifically instructed otherwise.
	"""
	velocity = 999
	if not auto:
		path_to_vel = os.path.join(malt.sd,
					   'malt90_velocities_year1.txt')
		f = open(path_to_vel,'r')
		for line in f:
			if line.split()[0].strip() == source:
				velocity = float(line.split()[1])

		f.close()
	if velocity == 999 or auto:
		velocity = identify_velocity(source,direction=direction)
	return(velocity)

def identify_velocity(source,minchan=200,maxchan=3896,sig=5,direction=None):
	"""Identify a source velocity based on HCO+
	Later I want to add HNC
	"""
	lines = ["hcop","hnc"]
	best_vel = {"hcop":0,"hnc":0}
	if direction:
		direction = direction+"_"
	for line in lines:
		try:
     			infile = get_filename(source,line,direction=direction)
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
				smoothspec = idl_stats.smooth(spec,
					     window_len=sf,window='hamming')
				mean,sigma = idl_stats.iterstat(smoothspec)
				goodsignal = np.where(smoothspec > 
						      threshold*sigma,1,0)
				goodsignal = scipy.ndimage.binary_erosion(
					     goodsignal,structure=np.ones(
					     channel_trim))
				maskedsignal = goodsignal*smoothspec
				max_chan[x,y] = np.argmax(maskedsignal)
		max_chan = np.extract(max_chan > 0,max_chan)
		best_chan = int(np.median(max_chan))+minchan
		best_vel[line] = vel[best_chan]/1000.
	ind_vels = []
	for line in lines:
		ind_vels.append(best_vel[line])
	ind_vels = np.array(ind_vels)
	return(np.median(ind_vels))

def do_source(source,lines,direction=None,auto=False):
	print("Sourcename: "+source)
	central_velocity = get_velocity(source,auto,direction)
	print("Center Velocity = "+str(central_velocity)+" km/s")
	create_basic_directories(lines)
	#create_output_directories(source,lines)
       	if direction:
	       	direction = direction+"_"
	for line in lines:
		infile = get_filename(source,line,direction)
		#print(infile)
		a = infile[:-9]
		out_base = a.replace("gridzilla","mommaps")
		out_dir = source+"_"+direction+line+"_mommaps"
		try:
			output_dir = os.path.join(malt.data_dir,
						 "mommaps",line,out_dir)
			os.mkdir(output_dir)
		except OSError:
			pass
		make_moment_maps(infile,out_base,output_dir,
				 central_velocity=central_velocity*1000)

def create_basic_directories(lines):
	"""Create subdirectories under mommaps"""
	for line in lines:
		try:
			moment_dir = os.path.join(malt.data_dir,"mommaps",line)
			os.mkdir(moment_dir)
		except OSError:
			pass
	
	
def get_filename(source,line,direction=None):
	"""Assumes direction has a _ at the end"""
	if not direction:
		filename = source+"_"+line+"_MEAN.fits"
	else:
		filename = source+"_"+direction+line+"_MEAN.fits"
	full_path = os.path.join(malt.data_dir,"gridzilla",line,filename)
	return(full_path)

def calculate_moments(d,minchan=False,maxchan=False,
		      vel=False,bestmask=False,mask=False):
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
	#print(minchan)
	#print(maxchan)
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
				mom = idl_stats.wt_moment(vel[ind],
						fullspec[ind],
						errors = np.zeros(ind.size)
							 +sigma)
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

def make_moment_maps(infile,out_base,output_dir,central_velocity=False,
		     second=False):
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
	

	#maps = do_predetermined_velocity(central_velocity,vel,hdout,n_edge,nchan,d,n_pad = 125)
	#save_maps(maps,hdout,out_base,out_dir,vel,minchan,maxchan,vwidth,"fullvel")
	
	maps = do_predetermined_velocity(central_velocity,vel,hdout,n_edge,
					 nchan,d,n_pad = 75)
	save_maps(maps,hdout,out_base,output_dir,vel,minchan,maxchan,vwidth,
		  "medvel")

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
	print("Velocity Integration: "+str(vmin)+' to '+str(vmax)+ 'km/s')
	hdout.update('VMIN',vmin,'KM/S')
	hdout.update('VMAX',vmax,'KM/S')
	mask = np.zeros(d.shape,dtype=np.int)
	mask[minchan:maxchan,...] = 1
	bestmask = mask[...,0,0]
	maps = calculate_moments(d,minchan,maxchan,vel,
				 bestmask=bestmask,mask=mask)
	return(maps)

def save_maps(maps,hdout,out_base,out_dir,vel,minchan,maxchan,vwidth,name_mod):
	
	(head,tail) = os.path.split(out_base)
	out_base2 = os.path.join(out_dir,tail)
#	out_temp = os.path.join(out_base2,tail)
	#print("Base for output")
	#print(out_base2)
#	print(out_temp)
	out_base = out_base2
	print("Saving maps to: "+out_base)
	#This trims out sources with sigma_v > 1000 km/s
	badind = np.where((maps['errmn'] > 1e6) | (maps['errsd'] > 1e6)) 
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
	
	ob = out_base
	maps['intint'] = maps['intint']*vwidth/1e3
	pyfits.writeto(ob+'mom1'+'.fits',maps['mean']/1e3,hdout,clobber=True)
	pyfits.writeto(ob+'mom2'+'.fits',maps['sd']/1e3,hdout,clobber=True)
	pyfits.writeto(ob+'err1'+'.fits',maps['errmn']/1e3,hdout,clobber=True)
	pyfits.writeto(ob+'err2'+'.fits',maps['errsd']/1e3,hdout,clobber=True)
	hdout.update('BUNIT','NONE')

	pyfits.writeto(ob+'npix'+'.fits',maps['npix'],clobber=True)
	pyfits.writeto(ob+'snr0'+'.fits',maps['intint']/(maps['error']*
		np.sqrt(maps['npix'])*vwidth/1e3),hdout,clobber=True)
	
	hdout.update('BUNIT','K.KM/S')
	
	pyfits.writeto(ob+'mom0'+'.fits',maps['intint'],hdout,clobber=True)
	pyfits.writeto(ob+'err0'+'.fits',maps['error']*np.sqrt(maps['npix'])*
					vwidth/1e3,hdout,clobber=True)
	pyfits.writeto(ob+'emap'+'.fits',maps['error']*
					vwidth/1e3,hdout,clobber=True)



def main():
	pass

if __name__ == '__main__':
	main()

