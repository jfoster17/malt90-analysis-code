#!/usr/bin/env python
#encoding: utf-8
"""
Examine a calibration file. This is a single-pointing of G301 (or G337)
used to verify setup and flux calibration. By default this just takes
the most recent calibration file of G301 for the current UT.

Options
-n : Night   -- specify which night (UT) to use when finding the most
                recent calibration file. Particularly useful if 
	        observations span the UT date change.
-c : CalFile -- specify a particular calibration file to use. This 
                is the raw UT date/time .rpf file.
-o : Object  -- specify which calibrator object to use. 
                Valid options: G301 (default) or G337
-a : All     -- Combine all calibration files for a given source
                to produce the highest SNR spectra. NOT WORKING
-h : Help    -- Display this help 
"""

import os,sys,getopt,glob,shutil,math,datetime
from subprocess import *
import numpy as np
import ReduceLog
import reduce_malt
import malt_params as malt
#This ugly stuff gets matplotlib 1.0.0 to be loaded
#sys.path[0] = '/usr/local/lib/python2.5/site-packages'
#sys.path[11] = ''
#sys.path[14] = ''
import matplotlib as mpl
#reload(mpl)
#...but it still does not work
#import matplotlib as mpl
mpl.use('Agg')
#print(mpl.__version__)
#import matplotlib.pylab as pylab

#import matplotl.pyplot as plt
#from asap import *
import asap
import pickle
reload(mpl)

def main():
	try:
		opts,args = getopt.getopt(sys.argv[1:],"n:c:o:ah")
	except getopt.GetoptError,err:
		print(str(err))
		print(__doc__)
		sys.exit(2)
	# Set default parameters #
	now_utc = datetime.datetime.utcnow()
        date   = datetime.date(now_utc.year,now_utc.month,now_utc.day).isoformat()
	do_all  = False
	filename = None
	source = "G301cal"
	ignore_date = True
	# Parse command line arguments #
	for o,a in opts:
		if o == "-n":
		       date = a
		       ignore_date = False
		elif o == "-c":
			filename = a
		elif o == "-o":
			object = a
		elif o == "-a":
			do_all = True
		elif o == "-h":
			print(__doc__)
			sys.exit(1)
		else:
			assert False, "unhandled option"
			print(__doc__)
			sys.exit(2)
	# Copy file over #
	redlog = ReduceLog.ReduceLog()
	redlog.update(date,in_middle_of_obs = True)
	if filename:
		id,source = redlog.get_name(filename)
	else:
		filename,source = redlog.find_latest_calibration_on_date(date,ignore_date=ignore_date)
	
	if not filename:
		print("No Calibration File Found")
		sys.exit(2)
        #Need to figure out how to handle "all" option well
	source = source.capitalize()
	source = source.rstrip("_")
	#Return the number of the cal file. If we have not done this source
	#Before it will return 0 and trigger the if
	already_done = redlog.check_cal(filename) 
	print(filename)
	print(already_done)
	#print(3+"bo")
#already_done = False
	if ignore_date:
		date = filename[0:10]
		print("Doing file from: "+date)
	if not already_done:
		current_files = glob.glob(malt.cal_dir+source+"*.rpf")
		index = 1
		for file in current_files:
			if date in file:
				index+=1
		renamed_file = source+"_"+date+"_"+str(index)+".rpf"
		shutil.copyfile(malt.source_dir+filename,malt.cal_dir+renamed_file)
       		redlog.mark_cal(filename,index) #Update redlog with index in rename slot
	else:
		renamed_file = source+"_"+date+"_"+str(already_done)+".rpf"

	if do_all and (source == "G301cal"):
		filelist = glob.glob(malt.cal_dir+"G301cal*.rpf")
		tempfile = filelist[0]
		curr = asap.scantable(tempfile,average=True)
		for filename in filelist:
			temp = asap.scantable(filename,average=True)
			curr = asap.merge(curr,temp)
		s = curr
		file = "G301cal_full"
	elif do_all and (source == "G337cal"):
		filelist = glob.glob(malt.cal_dir+"G337cal*.rpf")
		tempfile = filelist[0]
		curr = asap.scantable(tempfile,average=True)
		for filename in filelist:
			temp = asap.scantable(filename,average=True)
			curr = asap.merge(curr,temp)
       		s = curr
		file = "G337cal_full"
	else:
		s = asap.scantable(malt.cal_dir+renamed_file, average=True)
	cal_name = renamed_file.rstrip(".rpf")
	f_avt = plot_all_ifs(s,source,cal_name)
	dataline = fit_lines(f_avt,cal_name,source)
	update_database(source,dataline)
	plot_context(source,cal_name)
	print_report(source,cal_name)

def prep_scantable(s):
	# Set up basic ASAP stuff #
	lines,freqs,ifs = reduce_malt.setup_lines()
       	s.set_doppler('RADIO')
	s.set_freqframe('LSRK')
	q=s.auto_quotient()
	q.freq_align(insitu=True)
	q_pol=q.average_pol(weight='tsys')
	q_avt = q_pol.average_time(weight='tintsys',align=False)
	f_avt=q_avt.auto_poly_baseline(insitu=False,order=1,edge=[300])
	f_avt.smooth('gaussian',11,insitu=True)
	f_avt.set_restfreqs(freqs=freqs,unit="MHz")
	f_avt.set_unit("km/s")
	return(f_avt)

def plot_all_ifs(s,source,cal_name):
	"""Plot all 16 IFs for a given scantable.
	Return the averaged and set up scantable for other use.
	"""

	lines,freqs,ifs = reduce_malt.setup_lines()

	f_avt = prep_scantable(s)
	asap.plotter.set_legend(mode=-1)
	asap.plotter.set_colors('blue')
	asap.plotter.set_abcissa(fontsize=18)
	asap.plotter.set_ordinate(fontsize=18)
	asap.plotter.set_title(title=".",fontsize=4)
	asap.plotter.set_range(-120,200,-0.25,3.5) #G301cal
        #plotter.set_range(-120,200,-0.25,1.0) #Nessie

	calfolder = malt.cal_dir+"CalFolder"
	try:
		os.mkdir(calfolder)
	except OSError:
		pass
	#Plot each line and approximate central velocity
	for ifno,line in enumerate(lines):
		sel1 = asap.selector()
		sel1.set_ifs(ifno)
		f_avt.set_selection(sel1)
		asap.plotter.plot(f_avt)
		asap.plotter.text(100,1,line,fontsize=24)
		if (source == "G301cal"):
			asap.plotter.axvline(x=-42.7,ymin=0,ymax=1,color='r') #G301cal
		elif (source == "G337cal"):
			asap.plotter.axvline(x=-62.5,ymin=0,ymax=1,color='r') #G337cal
	        #plotter.axvline(x=-39,ymin=0,ymax=1,color='r') #Nessie
       		asap.plotter.save(calfolder+'/IF'+str(ifno).zfill(2)+'.eps')

	# Montage the individual images together #
	outname = malt.cal_dir+cal_name+".png"
       	p = Popen(["montage",calfolder+"/IF*.eps", "-tile", "4x4", \
				   "-geometry", "+0+0", outname])
	p.wait()
	shutil.copy(outname,malt.ver_dir)
	return(f_avt)

def fit_lines(f_avt,cal_name,source):
	lines,freqs,ifs = reduce_malt.setup_lines()

	dataline = np.zeros(1,dtype=[('name','a20'),('tsys','f4'),('elev','f4'),
				 ('n2hp','f4',6),('hnc','f4',6),
				 ('hcop','f4',6),('hcn','f4',6)])
	dataline[0]['name'] = cal_name
	fitlines = ["hcop","hnc","n2hp","hcn"]
	print(source)
	if (source.startswith("G301cal")):
		vbase = -42.7
	elif (source.startswith("G337cal")):
		vbase = -62.5
	for ifno,line in enumerate(lines):
		if line in fitlines:
			print(line)
			sel1 = asap.selector()
			sel1.set_ifs(ifno)
			f_avt.set_selection(sel1)
			dataline[0]['tsys'] = f_avt.get_tsys()[0]
			dataline[0]['elev'] = f_avt.get_elevation()[0]
			g = asap.fitter()
			g.set_scan(f_avt)
			if line == "hcop" or line == "hnc":
				n_gauss = 1
				g.set_function(gauss=1)
				g.set_gauss_parameters(2,vbase,5)
			elif line == "n2hp":
				n_gauss = 3
				g.set_function(gauss=3)
				g.set_gauss_parameters(1,vbase,4,component=0)
				g.set_gauss_parameters(0.5,vbase-9,4,component=1)
				g.set_gauss_parameters(0.5,vbase+6,4,component=2)
			elif line == "hcn":
				n_gauss = 3
				g.set_function(gauss=3)
				g.set_gauss_parameters(1,vbase,4,component=0)
				g.set_gauss_parameters(0.5,vbase-5,4,component=1)
				g.set_gauss_parameters(0.5,vbase+3,4,component=2)
			failed=False
			try:
				g.fit()
			except:
				failed=True
			if not failed:
				res = g.get_parameters()
			area_err = 1. #Dummy area error
			print(res['params'][0],res['errors'][0])
			print(res['params'][1],res['errors'][1])
			try:
				data_line = [res['params'][0],res['errors'][0],
					     res['params'][1],res['errors'][1],
					     g.get_area(),area_err]
				for j,entry in enumerate(data_line):
					dataline[0][line][j] = data_line[j]
			except IndexError:
				print("Failed Fit")
				failed = True
			if failed == False:
				asap.rcParams['plotter.gui'] = False
				fname = malt.cal_dir+cal_name+"_no_smooth_"+line+".png"
				g.plot(residual=True,filename=fname)
				shutil.copy(fname,malt.ver_dir)
	return(dataline)		

def read_database(source):
	cal_data_file = open(malt.cal_dir+source+"_caldbase.pkl",'rb')
	cal_data      = pickle.load(cal_data_file)
	cal_data_file.close()
	return(cal_data)

def write_database(source,cal_data):
       	cal_data_file = open(malt.cal_dir+source+"_caldbase.pkl",'wb')
	pickle.dump(cal_data,cal_data_file)
	cal_data_file.close()

def update_database(source,dataline):
	#Open current database
	cal_data = read_database(source)
        #Add dataline to database (check if it already exists)
	#Super inefficient way to do this
	length = len(cal_data)
	print(length)
	make_new = False
	if make_new:
		cal_data = dataline
		write_database(source,cal_data)
		return

	if dataline['name'] in cal_data['name']:
		for i in range(length):
			if cal_data[i]['name'] == dataline['name']:
				cal_data[i] = dataline
	else:
		new = np.zeros(length+1,dtype=[('name','a20'),
				 ('tsys','f4'),('elev','f4'),
				 ('n2hp','f4',6),('hnc','f4',6),
				 ('hcop','f4',6),('hcn','f4',6)])
		new[0:-1] = cal_data
		new[-1] = dataline
		cal_data = new
#	print(cal_data)
	#Save out database.
	write_database(source,cal_data)

def plot_context(source,cal_name):
	"""Plot the latest addition to the databse in relation to all others"""
	temp_cal_data = read_database(source)
	cal_data = np.sort(temp_cal_data, order='name')
#	print(cal_data)
	new_data = cal_data[cal_data['name'] == cal_name]
	old_data = cal_data[cal_data['name'] != cal_name]
	n_old = len(old_data)
#	print(new_data['hnc'][0][0])
       	fitlines = ["hcop","hnc","n2hp","hcn"]
	colors = ["red","green","blue","purple"]
	nums = [4,3,2,1]
	for line,ccolor,i in zip(fitlines,colors,nums):
		if i == 4:
			ax1 = mpl.pylab.subplot(4,1,i)
			mpl.pylab.ylabel("Peak T [K]",fontsize=9)
		else:
			ax2 = mpl.pylab.subplot(4,1,i,sharex=ax1)
			mpl.pylab.setp(ax2.get_xticklabels(),visible=False)
		dates2 = old_data['name']
		dates = [date[-12:] for date in dates2]
		old_line = old_data[line]
		prior_mean = np.median(old_line[:,0])
		prior_range = np.std(old_line[:,0])
		yup = prior_mean+prior_range
		ydo = prior_mean-prior_range
		#print sys.path
		mpl.pylab.fill([-1,-1,n_old+1,n_old+1],[ydo,yup,yup,ydo],fc=ccolor,alpha=0.1)
      		yup = prior_mean+prior_range*2
		ydo = prior_mean-prior_range*2
		mpl.pylab.fill([-1,-1,n_old+1,n_old+1],[ydo,yup,yup,ydo],fc=ccolor,alpha=0.1)
		mpl.pylab.plot(old_line[:,0],'o',color=ccolor,ms=4)
		mpl.pylab.plot([n_old],[new_data[line][0][0]],'+',mec=ccolor,mfc=ccolor,ms=10)
		mpl.pylab.title(line,fontsize=11)
		mpl.pylab.xlim(-1,n_old+1) 
		locs,lables = mpl.pylab.xticks()
		xpos = np.arange(-1,n_old+1)
		xlab = ['']+list(dates)+[new_data['name'][0][-12:]]+['']
		if len(xpos) % 2 == 0:
			st = 1
		else:
			st = 0
		mpl.pylab.xticks(xpos[st::2],xlab[st::2],rotation=30,fontsize=5,ha='right')
#		locs,labels = pylab.yticks()
		mpl.pylab.yticks(fontsize=9)
	mpl.pylab.savefig(malt.ver_dir+source+"LatestCal.png")
	mpl.pylab.savefig(malt.ver_dir+"LatestCal.png")


def print_report(source,cal_name):
	"""Print out a text report about the quality of this calibration file"""
	cal_data = read_database(source)
#	print(cal_data)
	new_data = cal_data[cal_data['name'] == cal_name]
	old_data = cal_data[cal_data['name'] != cal_name]
       	fitlines = ["hcop","hnc","n2hp","hcn"]
	histpeak = {'hcop':'3.0 +/- 0.5 K','hcn':'2.2 +/- 0.3','hnc':'1.8 +/- 0.4','n2hp':'1.2 +/- 0.8'}
	for line in fitlines:
		old_line = old_data[line]
		median = np.median(old_line[:,0])
		stan =  np.std(old_line[:,0])
		histpeak[line] = "%2.1f +/- %2.1f" % (median,stan)

 	print("#################### Calibration file summary ###################")
	print("Filename -- "+cal_name+'.rpf')
	print("Look at /DATA/MALT_1/MALT90/data/cal/"+cal_name+'.png'+" to see the spectra")
	for line in fitlines:
		print("Fit paramters for: "+line)
		print("    Peak     = "+str(new_data[0][line][0]))
		print("    Normal Range is "+histpeak[line])
#		print("Velocity = "+str(data[0][line][2]))
	print("#################### Calibration file summary ###################")


if __name__ == '__main__':
	main()
