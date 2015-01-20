include 'gridzilla.g'

name := argv[2]

# Input directory.
read_dir := spaste('/DATA/MALT_1/MALT90/data/byhand/livedata/',name)

# Output directory.
write_dir := spaste('/DATA/MALT_1/MALT90/data/byhand/gridzilla/',name)


doif := [F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F]

infile := argv[6:len(argv)]
outfile := argv[5]
rfreq := as_float(argv[3])
ifno  := as_integer(argv[4])+1
#print argv
#print ifno
doif[ifno] := T
#print rfreq

# Velocity/frequency information
vel1 := -1000
vel2 := 1000

# Start the gridder.
lsg := gridzilla(remote = T,
                 autosize = T,
                 pixel_width = 0.15,
                 pixel_height = 0.15,
                 rangeSpec = 'velocity',
                 startSpec = vel1,
                 endSpec = vel2,
                 restFreq = rfreq,
                 IFsel = doif,
                 pol_op = 'A&B',
                 spectral = T,
                 continuum = F,
                 baseline = F,
                 projection = 'SIN',
                 coordSys = 'galactic',
                 statistic = 'mean',
                 clip_fraction = 0.0,
                 tsys_weight = T,
                 beam_weight = 0,
                 beam_FWHM = 0.6,
                 beam_normal = F,
                 kernel_type = 'gaussian',
                 kernel_FWHM = 0.2,
                 cutoff_radius = 0.2,
                 storage = 33,
                 chan_err = 0.2,
                 directories = read_dir,
                 files = infile,
                 selection = ind(infile),
                 write_dir = write_dir,
                 p_FITSfilename=outfile,
                 spectype = 'VELO-XXX',
                 short_int = F)

# Process the data
lsg->go()
await lsg->finished
print 'Finished gridzilla processing.'

exit


