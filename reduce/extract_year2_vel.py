import os

f = open('MALT90Y2catalog.dat','r')
g = open('temp.txt','w')
for line in f:
    if not line.startswith("#"):
        source,type,lon,lat,vbest,vrms,vcode,kdnear,kdfar,derrp,derrm,probn,probf,decod,comment = line.split()
        print >>g,source,vbest
f.close()
g.close()
