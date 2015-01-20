import pyfits
import os
import numpy as np


def main():
    dir = "/DATA/MALT_1/MALT90/data/year1/sources/"
    year1sources = os.listdir(dir)
    print("Year1")
    find_noisy(year1sources,dir)
    
    print("Year2")
    dir = "/DATA/MALT_1/MALT90/data/year2/sources/"
    year2sources = os.listdir(dir)
    find_noisy(year2sources,dir)

    print("Year3")
    dir = "/DATA/MALT_1/MALT90/data/sources/"
    year3sources = os.listdir(dir)
    find_noisy(year3sources,dir)

def find_noisy(sources,dir):
    for source in sources:
#        print(source)
        try:
            longitude = float(source[1:4])
        except ValueError:
            longitude = 0.
        try:
            d,h = pyfits.getdata(os.path.join(dir,source,source+"_13c34s_mommaps",source+"_13c34s_emap.fits"),header=True)
            noise_13cs = np.median(np.nan_to_num(d.flatten()))
        except IOError:
            print("Failed to find "+source)
            noise_13cs = 0

        
        try:
            d,h = pyfits.getdata(os.path.join(dir,source,source+"_h13cop_mommaps",source+"_h13cop_emap.fits"),header=True)
            noise_h13cop = np.median(np.nan_to_num(d.flatten()))
        except IOError:
            noise_h13cop = 0
        print(source+":   "+str(noise_13cs)+"  "+str(noise_h13cop))

main()
