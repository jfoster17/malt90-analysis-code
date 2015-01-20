#!/usr/bin/env python
# encoding: utf-8
import malt_params as malt
import sys
import os
import Malt90SourceBasic
from sets import Set

def main():
    all_ids = get_all_ids(year_filter = "2011")
    all_sources = []
    finished_sources = []
    half_sources = []
    #print(all_ids)
    #print(len(all_ids))
    for source_id in all_ids:
        single_files = []
        list_files = []
        source_names = []
        source = Malt90SourceBasic.Malt90SourceBasic(source_id)
        all_sources.append(source)
        for filename in source.source_names:
            if "B" not in filename:
                single_files.append(filename)
                set_files = Set(single_files)
                list_files = list(set_files)
        if len(list_files) > 1:
            source.trimmed_files = list_files
            finished_sources.append(source)
        if len(list_files) == 1:
            source.trimmed_files = list_files
            half_sources.append(source)

    #print(finished_sources)
    #return(0)

    for source in finished_sources:
        #print(source.id,source.trimmed_files)
        name = source.trimmed_files[0].split('_')[0]
        source_names.append(name)
        #source.get_noise()
        #one_bad = False
        #for i,noise in enumerate(source.maps_noise):
        #    if noise > 0.9:
        #        one_bad = True
        #if one_bad:
        #    print(source.basic_name)
        #    print(source.maps_noise)
        #    print(source.trimmed_files)

    for source_name in source_names:
        print(source_name)
#    print(len(finished_sources))
 #   print(len(half_sources))
  #  for source in half_sources:
   #     print(source.id,source.trimmed_files)
    
            
def get_all_ids(year_filter = None):
    f = open(malt.log_location,'r')
    lines = f.readlines()
    all_ids = []
    for line in lines:
        if line.startswith("#"):
            pass
        else:
            en = line.split()
            try:
                raw_filename = en[0]
                atlasgal_id  = en[1]
                sourcename   = en[2]
                if (year_filter in raw_filename) and (int(atlasgal_id) != 0):
                    all_ids.append(atlasgal_id)
            except IndexError: #Deal with blank line
                print("Blank or corrupt line in reduction log")
    set_ids = Set(all_ids)
    list_ids = list(set_ids)
    f.close()
    return(list_ids)

if __name__ == '__main__':
    main()
