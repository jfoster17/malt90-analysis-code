import sys
import os
import Malt90SourceSNR

def main():
    output_file = "year1to3_snr_stats.txt"
    ff = open('year1to3_sourcelist.txt','r')
    gg = open(output_file,'w')
    print >>gg,"#Sourcename     n2hp_max    n2hp_num   hnc_max   hnc_num    hcop_max   hcop_num   hcn_max   hcn_num"
    for line in ff:
        name = line.strip()
        a = Malt90SourceSNR.Malt90SourceSNR(name)
        gg.write(name)
        a.get_all_snr()
        for line in ["n2hp","hnc","hcop","hcn"]:
            gg.write("    "+str(a.max_snr[line])+"    "+str(a.num_snr[line]))
        gg.write("\n")
    gg.close()
    ff.close()



if __name__ == '__main__':
    main()
