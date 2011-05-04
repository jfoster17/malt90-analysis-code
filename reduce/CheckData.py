#!/usr/bin/env python
# encoding: utf-8
"""
untitled.py

Created by Jonathan Foster on 2010-09-17.
Copyright (c) 2010 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import pyfits
import numpy as np

def main():
	examine_folder("/DATA/MALT_1/MALT90/data/old_moment_maps/hcop/")

def examine_folder(path):
	files = os.listdir(path)
	for filename in files:
		basename = filename.replace("_mom_maps",'')
		data = pyfits.getdata(path+filename+'/'+basename+".emap.fits")
		print(basename)
		#print(data)
		print(np.median(data,axis=None))

if __name__ == '__main__':
	main()

