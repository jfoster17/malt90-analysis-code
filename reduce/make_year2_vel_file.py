f = open('year2_sources.txt','r')
g = open('auto_vel.txt','r')

sources = []
for line in f:
    source = line.strip()
    sources.append(source)

vel_sources = []
vels = []
for line in g:
#    print(line)
    source,vel = line.split()
    vel_sources.append(source)
    vels.append(vel)

f.close()
g.close()

out_vels = []
for source in sources:
    found_flag = False
    for i,vel_source in enumerate(vel_sources):
        if source == vel_source:
            out_vels.append(vels[i])
            found_flag = True
            break
    if not found_flag:
        out_vels.append(0)

for source,vel in zip(sources,out_vels):
    print(source+"\t"+str(vel))

#print(sources)
#print(out_vels)
#print(len(sources))
#print(len(out_vels))
