#!/usr/X11R6/bin/python

from asap import *
import os,sys,getopt
from subprocess import *
import glob
import malt_params as malt
import ReduceLog
import shutil
import datetime
import numpy as np
import math

#Malt90 line parameters
lines = ["n2hp","13cs","h41a","ch3cn",\
	"hc3n","13c34s","hnc","hc13ccn",\
	"hcop","hcn","hnco413","hnco404",\
	"c2h","hn13c","sio","h13cop"]
freqs = [93173.772, 92494.303, 92034.475, 91985.316, \
         90978.989, 90926.036, 90663.572, 90593.059, \
         89188.526, 88631.847, 88239.027, 87925.238, \
         87316.925, 87090.850, 86847.010, 86754.330]

def main():
	try:
		opts,args = getopt.getopt(sys.argv[1:],"n:c:a")
	except getopt.GetoptError,err:
		print(str(err))
		print(__doc__)
		sys.exit(2)
	now_utc = datetime.datetime.utcnow()
        date   = datetime.date(now_utc.year,now_utc.month,now_utc.day).isoformat()
	do_all  = False
	filename = None
	for o,a in opts:
		if o == "-n":
		       date = a
		if o == "-c":
			filename = a
		if o == "-a":
			source = a
			do_all = True

	redlog = ReduceLog.ReduceLog()
	redlog.update(date,in_middle_of_obs = True)
	if filename:
		id,source = redlog.get_name(filename)
	else:
		filename,source = redlog.find_latest_calibration_on_date(date)
	#Need to figure out how to handle "all" option well

	source = source.capitalize()
	current_files = glob.glob(malt.cal_dir+source+"*.rpf")
	index = 1
	for file in current_files:
		if date in file:
			index+=1
	renamed_file = source+"_"+date+"_"+str(index)+".rpf"
	shutil.copyfile(malt.source_dir+filename,malt.cal_dir+renamed_file)
	#Copy file over


	if do_all and (source == "G301cal"):
		filelist = glob.glob(malt.cal_dir+"G301cal*.rpf")
		tempfile = filelist[0]
		curr = scantable(tempfile,average=True)
		for filename in filelist:
			temp = scantable(filename,average=True)
			curr = merge(curr,temp)
		s = curr
		file = "G301cal_full"
	elif do_all and (source == "G337cal"):
		filelist = glob.glob(malt.cal_dir+"G337cal*.rpf")
		tempfile = filelist[0]
		curr = scantable(tempfile,average=True)
		for filename in filelist:
			temp = scantable(filename,average=True)
			curr = merge(curr,temp)
       		s = curr
		file = "G337cal_full"
	else:
		s = scantable(malt.cal_dir+renamed_file, average=True)

#s.save(name="AllScans.rpf",format='ASAP')

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
	plotter.set_legend(mode=-1)
	plotter.set_range(-120,200,-0.25,3.5) #G301cal
#plotter.set_range(-120,200,-0.25,1.0) #Nessie

	calfolder = malt.cal_dir+"CalFolder"
	try:
		os.mkdir(calfolder)
	except OSError:
		pass

	plotter.set_colors('blue')
	plotter.set_abcissa(fontsize=18)
	plotter.set_ordinate(fontsize=18)
	plotter.set_title(title=".",fontsize=4)

	for ifno,line in enumerate(lines):
		sel1 = selector()
		sel1.set_ifs(ifno)
		f_avt.set_selection(sel1)
		plotter.plot(f_avt)
		plotter.text(100,1,line,fontsize=24)
		if (source == "G301cal"):
			plotter.axvline(x=-42.7,ymin=0,ymax=1,color='r') #G301cal
		elif (source == "G337cal"):
			plotter.axvline(x=-62.5,ymin=0,ymax=1,color='r') #G337cal
	#plotter.axvline(x=-39,ymin=0,ymax=1,color='r') #Nessie
       		plotter.save(calfolder+'/IF'+str(ifno).zfill(2)+'.eps')

	outname = malt.cal_dir+renamed_file.rstrip(".rpf")+".pdf"
       	p = Popen(["montage",calfolder+"/IF*.eps", "-tile", "4x4", \
				   "-geometry", "+0+0", outname])
	p.wait()


	data = np.zeros(1,dtype=[('name','a20'),('tsys','f4'),('elev','f4'),('n2hp','f4',6),('hnc','f4',6),('hcop','f4',6),('hcn','f4',6)])

	fitlines = ["hcop","hnc","n2hp","hcn"]
#Plot each line and approximate central velocity
	for ifno,line in enumerate(lines):
		if line in fitlines:
			sel1 = selector()
			sel1.set_ifs(ifno)
			f_avt.set_selection(sel1)
			g = fitter()
			g.set_scan(f_avt)
			if line == "hcop" or line == "hnc":
				n_gauss = 1
				g.set_function(gauss=1)
				g.set_gauss_parameters(2,-50,5)
			elif line == "n2hp":
				n_gauss = 3
				g.set_function(gauss=3)
				g.set_gauss_parameters(1,-41,4,component=0)
				g.set_gauss_parameters(0.5,-50,4,component=1)
				g.set_gauss_parameters(0.5,-35,4,component=2)
			elif line == "hcn":
				n_gauss = 3
				g.set_function(gauss=3)
				g.set_gauss_parameters(1,-41,4,component=0)
				g.set_gauss_parameters(0.5,-48,4,component=1)
				g.set_gauss_parameters(0.5,-37,4,component=2)

			g.fit()
			failed = False
			res = g.get_parameters()
			chi2 = g.get_chi2()
			if chi2 > 1:
				errfac = 1.
	#errfac = math.sqrt(chi2)
			else:
				errfac = 1.

			try:
				if n_gauss == 3:
					area_err = math.sqrt((res['errors'][0]*errfac/res['params'][0])**2+(res['errors'][2]*errfac/res['params'][2])**2
							     +(res['errors'][3]*errfac/res['params'][3])**2+(res['errors'][5]*errfac/res['params'][5])**2
							     +(res['errors'][6]*errfac/res['params'][6])**2+(res['errors'][8]*errfac/res['params'][8])**2)
				else:
					area_err = math.sqrt((res['errors'][0]*errfac/res['params'][0])**2+(res['errors'][2]*errfac/res['params'][2])**2)
				data_line = [res['params'][0],res['errors'][0]*errfac,res['params'][1],res['errors'][1]*errfac,g.get_area(),area_err]
	
       				for j,entry in enumerate(data_line):
					data[0][line][j] = data_line[j]
			except IndexError:
				print("Failed Fit")
				failed = True
			if failed == False:
				rcParams['plotter.gui'] = False
				g.plot(residual=True,filename=malt.cal_dir+renamed_file.strip('.rpf')+"_no_smooth_"+line+".png")
	histpeak = {'hcop':'3.0 +/- 0.5 K','hcn':'2.2 +/- 0.3','hnc':'1.8 +/- 0.4','n2hp':'1.2 +/- 0.8'}
 	print("#################### Calibration file summary ###################")
	print("Filename -- "+renamed_file)
	print("Look at /DATA/MALT_1/MALT90/data/cal/"+renamed_file.replace('.rpf','.pdf')+" to see the spectra")
	for line in fitlines:
		print("Fit paramters for: "+line)
		print("    Peak     = "+str(data[0][line][0]))
		print("    Normal Range is "+histpeak[line])
#		print("Velocity = "+str(data[0][line][2]))
	print("#################### Calibration file summary ###################")

#Install gqview
#	p = Popen(["xpdf",outname])
	#p.wait()

if __name__ == '__main__':
	main()
