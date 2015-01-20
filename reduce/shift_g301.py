import pyfits
import os,sys
import reduce_malt
import numpy as np

lines,freqs,ifs = reduce_malt.setup_lines()

for line in lines:
    for i in range(4): #Ignore/do not modify the last two
        fname = '/DATA/MALT_1/MALT90/data/byhand/livedata/'+line+'/G300.969+01.145_'+str(i)+'_'+line+'.sdfits'
        d,h = pyfits.getdata(fname,header=True)
        if i in (0,1): #2011-05-06
            shiftx = 0.00243
            shifty = -0.001199
        elif i in (2,3): #2011-08-22
            shiftx = 0.0066
            shifty = 0.00039
        d.field('CRVAL3')[:]  = d.field('CRVAL3')[:]+shiftx
        d.field('OBJ-RA')[:]  = d.field('OBJ-RA')[:]+shiftx
        d.field('CRVAL4')[:]  = d.field('CRVAL4')[:]+shifty
        d.field('OBJ-DEC')[:] = d.field('OBJ-DEC')[:]+shifty
        fname = '/DATA/MALT_1/MALT90/data/byhand/livedata/'+line+'/G300.969+01.145_'+str(i)+'_'+line+'.sdfits'
        pyfits.writeto(fname,d,h,clobber=True)
        
        
        
        
