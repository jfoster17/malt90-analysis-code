#!/usr/X11R6/bin/python

from asap import *
import os,sys
from subprocess import *
import glob
import malt_params as malt

malt.cal_dir

#Malt90 line parameters
lines = ["n2hp","13cs","h41a","ch3cn",\
	"hc3n","13c34s","hnc","hc13ccn",\
	"hcop","hcn","hnco413","hnco404",\
	"c2h","hn13c","sio","h13cop"]
freqs = [93173.772, 92494.303, 92034.475, 91985.316, \
         90978.989, 90926.036, 90663.572, 90593.059, \
         89188.526, 88631.847, 88239.027, 87925.238, \
         87316.925, 87090.850, 86847.010, 86754.330]


file = sys.argv[1]

if file == "all":
	filelist = glob.glob(data_dir+"G301cal*.rpf")
	tempfile = filelist[0]
	curr = scantable(tempfile,average=True)
	for filename in filelist:
		temp = scantable(filename,average=True)
		curr = merge(curr,temp)
	s = curr
	file = "G301cal_full"
elif file == "all337":
	filelist = glob.glob(data_dir+"G337cal*.rpf")
	tempfile = filelist[0]
	curr = scantable(tempfile,average=True)
	for filename in filelist:
		temp = scantable(filename,average=True)
		curr = merge(curr,temp)
	s = curr
	file = "G337cal_full"
else:
	s = scantable(data_dir+file, average=True)

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

calfolder = data_dir+"CalFolder"
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
	if file.startswith("G301cal"):
		plotter.axvline(x=-42.7,ymin=0,ymax=1,color='r') #G301cal
	elif file.startswith("G337cal"):
		plotter.axvline(x=-62.5,ymin=0,ymax=1,color='r') #G337cal
	#plotter.axvline(x=-39,ymin=0,ymax=1,color='r') #Nessie
	plotter.save(calfolder+'/IF'+str(ifno).zfill(2)+'.eps')

outname = data_dir+file.rstrip(".rpf")+".pdf"
p = Popen(["montage",calfolder+"/IF*.eps", "-tile", "4x4", \
	"-geometry", "+0+0", outname])
p.wait()

#Install gqview
p = Popen(["xpdf",outname])
p.wait()
