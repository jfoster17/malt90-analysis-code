#!/usr/bin/env python
# encoding: utf-8
"""
This program is supposed to rename the files that we incorrectly named in the original MALT90 year1 and year2 data releases.



"""

import sys
import os
import string
import shutil

lines = ["13c34s","c2h","h13cop","hc13ccn","hcn","hn13c","hnco404","n2hp","13cs","ch3cn","h41a","hc3n","hcop","hnc","hnco413","sio"]

#lines = ["13c34s","sio"]

def main():
    rename(old_name = "G314.194+00.269", new_name = "G302.032-00.061")
    rename(old_name = "G311.627+00.290", new_name = "G314.194+00.269")
    #rename(old_name = "G302.032-00.061", new_name = "G303.930-00.688")
    
def rename(old_name=None,new_name=None):
    
    #Walk through the gridzilla folder renaming all the FITS cubes.
    #This is easy and straightforward
    bd = "/DATA/MALT_1/MALT90/data/yearx/gridzilla/"
    for line in lines:
        oldfile = bd+line+"/"+old_name+"_"+line+"_MEAN.fits"
        print(">>"+oldfile)
        newfile = bd+line+"/"+new_name+"_"+line+"_MEAN.fits"
        print("<<"+newfile)
        os.rename(oldfile,newfile)
        
    #Now deal with the sources/ folder.
    #Remove the pre-existing sources/ folder if necessary
    #Move over the full directory
    #Step through the subdirectories, renaming files.
    #Don't know how to do the symlink though.
    bd = "/DATA/MALT_1/MALT90/data/yearx/sources/"
    try:
        shutil.rmtree(bd+new_name)
    except OSError:
        pass
    print(">>"+bd+old_name)
    print("<<"+bd+new_name)
    shutil.move(bd+old_name,bd+new_name)
    top_dir = bd+new_name    
    allfiles = os.listdir(top_dir)
    print(allfiles)
    for ifile in allfiles:
        if os.path.isdir(top_dir+"/"+ifile) == True:
            subfiles = os.listdir(top_dir+"/"+ifile)
            for jfile in subfiles:
                print(">>"+top_dir+"/"+ifile+"/"+jfile)
                print("<<"+top_dir+"/"+ifile+"/"+string.replace(jfile,old_name,new_name))
                os.rename(top_dir+"/"+ifile+"/"+jfile,top_dir+"/"+ifile+"/"+string.replace(jfile,old_name,new_name))
            print(">>"+top_dir+"/"+ifile)
            print("<<"+top_dir+"/"+string.replace(ifile,old_name,new_name))  
            shutil.move(top_dir+"/"+ifile,top_dir+"/"+string.replace(ifile,old_name,new_name))
        else:
            print(">>"+top_dir+"/"+ifile)
            print("<<"+top_dir+"/"+string.replace(ifile,old_name,new_name))
            if "MEAN" in ifile:
                try:
                    os.unlink(top_dir+"/"+ifile)
                except:
                    pass
                try:
                    line_name = ifile.split('_')[1]
                    print(line_name)
                    os.symlink("../../gridzilla/"+line_name+"/"+string.replace(ifile,old_name,new_name),top_dir+"/"+string.replace(ifile,old_name,new_name))
                except:
                    pass
            else:
                os.rename(top_dir+"/"+ifile,top_dir+"/"+string.replace(ifile,old_name,new_name))
            
            #This works for the masks
            #But it probably doesn't work for symlinks
    
if __name__ == '__main__':
    main()

