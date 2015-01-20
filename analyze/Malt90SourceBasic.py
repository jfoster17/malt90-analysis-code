
import malt_params as malt
import pyfits
import os
import numpy as np

class Malt90SourceBasic:
    """A class to store basic information about a Malt90Source.
    Can run on Draco, so limited imports required.
    Stores information about
    (1) When source was observed
    (2) What the noise properties are like
    (3) Which files went into the reduction. 
    """
    def __init__(self,atlasgal_id):
        self.id = atlasgal_id
        self.comments = self.get_comments()
        self.get_log_info()
        self.set_basic_name()
        
        self.data_dir = "/DATA/MALT_1/MALT90/data/sources/"+self.basic_name+"/"
        self.data_dirs = [self.data_dir+"year1",self.data_dir+"year2",self.data_dir+"year3"]


    def set_basic_name(self):
        self.basic_name = self.source_names[0][0:15]
        #print(self.basic_name)

    def get_log_info(self):
        self.location = malt.log_location
        f = open(self.location,'r')
        lines = f.readlines()
        self.raw_filenames = []
        self.source_names  = []
        for line in lines:
            if line.startswith("#"):
                pass
            else:
                en = line.split()
                try:
                    raw_filename = en[0]
                    atlasgal_id  = en[1]
                    sourcename   = en[2]
                    if atlasgal_id == self.id:
                        self.raw_filenames.append(raw_filename)
                        self.source_names.append(sourcename)
                except IndexError: #Deal with blank line
                    print("Blank or corrupt line in reduction log")
        f.close()

    def get_comments(self):
        """Get freeform comments about a source from text database."""
        pass
    
    def robust_get_moment_map(self,molecule,moment):
        """
        Terrible kludge to get moment maps out of all three years
        """
        try:
            path = self.get_moment_map_path(self.data_dirs[2],molecule,moment)
            data = np.nan_to_num(pyfits.getdata(path))
        except IOError:
            try: 
                path = self.get_moment_map_path(self.data_dirs[1],molecule,moment)
                data = np.nan_to_num(pyfits.getdata(path))
            except IOError:
                try: 
                    path = self.get_moment_map_path(self.data_dirs[0],molecule,moment)
                    data = np.nan_to_num(pyfits.getdata(path))
                except IOError:
                    print("No such moment map found!")
                    return(0)
        return(data)

    def get_moment_map_path(self,dir,molecule,moment):
        mol = "_"+molecule+"_"
        a = os.path.join(dir,self.name+mol+"mommaps",self.name+mol+moment+".fits")
        print(a)
        return(a)

    def get_all_snr(self):

        self.max_snr = {"n2hp":0,"hnc":0,"hcop":0,"hcn":0}
        self.num_snr = {"n2hp":0,"hnc":0,"hcop":0,"hcn":0}
        for line in ['n2hp','hnc','hcop','hcn']:
            self.max_snr[line] = self.get_snr(line)
            self.num_snr[line] = self.get_num_snr_in_center(line,threshold=5)
            
    def get_snr(self,line):
        """
        Get the maximum SNR in the map. 
        """
        snr = self.robust_get_moment_map(line,'snr0')
        max_val = np.max(snr[8:18,8:18])
        return(max_val)

    def get_num_snr_in_center(self,line,threshold=5):
        snr = self.robust_get_moment_map(line,'snr0')
        above_thresh = np.ma.masked_less(snr[8:18,8:18],threshold)
        num = np.ma.count(above_thresh)
        return(num)

    def get_noise(self):
        """Estimate the noise for a source (both GLon and GLat) from 13c34s cube."""
        all_maps = self.trimmed_files
#        both_maps = [self.basic_name+"_GLat_13c34s_MEAN.fits",self.basic_name+"_GLon_13c34s_MEAN.fits"]
        self.maps_noise = []
        for i,one_map in enumerate(all_maps):
            try:
                d,h = pyfits.getdata(os.path.join(malt.data_dir,'gridzilla','13c34s',one_map+"_13c34s_MEAN.fits"),header=True)
                cen_spec = d[200:3896,13,13]
                noise_estimate = np.std(cen_spec)
                self.maps_noise.append(noise_estimate)
            except IOError:
                print("Failed to find "+one_map)
                self.maps_noise.append(999.)


            #print(map_noise)
        #self.GLat_noise = GLat_temp
        #self.GLon_noise = GLon_temp
