# Input files given on the command line (wildcards supported).
#
# Usage: glish -l ldata_if3.g -plain <infile1> <infile2>...
# Usage: glish -l ldata_if_16.g -plain <infile1> <infile2>...

include 'livedatareducer.g'

#Normal line setup
lines := ["n2hp","13cs","h41a","ch3cn","hc3n","13c34s","hnc","hc13ccn","hcop","hcn","hnco413","hnco404","c2h","hn13c","sio","h13cop"]
#First night/backward IFS
#lines := ["ch3cn","h41a","13cs","n2hp","hc13ccn","hnc","13c34s","hc3n","hnco404","hnco413","hcn","hcop","c2h","hn13c","sio","h13cop"]

# Input directory.
read_dir := '/DATA/MALT_1/MALT90/data/renamed'
files := argv[3:len(argv)]

# Output directory.
write_dir := '/DATA/MALT_1/MALT90/data/livedata'

# Start and stop channels
startc := 1
endc :=  4096

# Instantiate a livedata reducer.
ldred := reducer(bandpass   = T,
                 monitor    = F,
                 writer     = T,
                 stats      = F,
                 read_dir   = read_dir,
                 write_dir  = write_dir)

# Create a GUI so we can see what's going on.
ldred->showgui()
t := client("timer -oneshot", 3.0)
await t->ready

# Configure the bandpass client.
ldred.bandpass->setparm([smoothing='HANNING', 
                 method='REFERENCED',
                 xbeam=F,
                 fit_order=2,
                 doppler_frame='LSRK',
                 rescale_axis=T,
                 continuum=F,
                 chan_mask = [1,250,0,0,0,0,0,0,3896,4096],
                 estimator='MEAN'])


# Process a batch of files through livedata.
for (read_file in files) {
  print 'Begin processing' 
  if (read_file !~ m/\.rpf/) next
  print 'Processing', read_file
  doif := [F,F,F,F,F,F,F,F,F,F,F,F,F,F,F,F]
  ifno := 9
  print 'Working on IF', ifno
  doif[ifno] := T
  rootname := read_file ~ s/(\.rpf|\.rpf)//
  write_file := spaste(lines[ifno],'/',rootname,'_',lines[ifno])
  ldred.reader->setparm(IFsel=doif)
  ldred.reader->setparm(startChan=startc)
  ldred.reader->setparm(endChan=endc)
    #directory := spaste('sdfits/',lines[ifno])
    #ldred.writer->setparm(write_dir = )
  ldred->start([read_file=read_file,write_file=write_file])
  await ldred->finished
  doif[ifno] := F
}

# Finished with livedata.
print 'Closing down livedata...'
ldred->terminate()
await ldred->done
print 'Finished livedata processing.'

exit

