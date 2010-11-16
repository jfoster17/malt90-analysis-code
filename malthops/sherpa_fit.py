#!/usr/bin/env python 
import sys 
sys.path.append('/usr/bin') 
sys.path.append('/usr/lib/python25.zip') 
sys.path.append('/usr/lib/python2.5') 
sys.path.append('/usr/lib/python2.5/plat-linux2') 
sys.path.append('/usr/lib/python2.5/lib-tk') 
sys.path.append('/usr/lib/python2.5/lib-dynload') 
sys.path.append('/usr/local/lib/python2.5/site-packages') 
sys.path.append('/usr/lib/python2.5/site-packages') 
sys.path.append('/usr/lib/python2.5/site-packages/Numeric') 
sys.path.append('/usr/lib/python2.5/site-packages/gst-0.10') 
sys.path.append('/var/lib/python-support/python2.5') 
sys.path.append('/usr/lib/python2.5/site-packages/gtk-2.0') 
sys.path.append('/var/lib/python-support/python2.5/gtk-2.0') 
sys.path.append('/var/lib/python-support/python2.5/IPython/Extensions') 
import pyfits as pf
import numpy as np
import qfit_tl_hops as qfit
res = qfit.fitfile("/mako3/MALT_1/MALT90/data/gridzilla/hnc/_hnc_MEAN.fits", ngauss=1, verbose=True,linkfwhm=False,linkpos=False,offsets=[0],maskimg=None,plotname=None) 
qfit.savefit(res,"_hnc_MEAN_sherpa") 
