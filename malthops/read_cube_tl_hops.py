#
# experimentation at reading in a FITS data cube using crates
# and getting the necessary metadata through hacky and horrible
# means.
#

#import pycrates
from pycrates import *
import subprocess
import numpy

def get_example_data():
    """Return x,y values from spectrum.dat.
    
    This doesn't really fit in read_cube.py but leave for now.
    """
    cr = read_file ("spectrum.dat[opt skip=7]")
    x = get_colvals(cr,0)
    y = get_colvals(cr,1)
    # we copy the values to avoid memory-management issues
    return (x.copy(), y.copy())

def read_cube(filename):
    """Read in a data cube and return the data, plus
    the velocity coordinates (i.e. the location of the Z
    pixels).

    This is complicated due to bugs and limitations in
    Crates in CIAO 4.0.
    """
    f = open('junk.txt', 'a')
    print >>f,"Starting Read"
    print >>f,filename
    try:
        cr = read_file(filename)
    except OSError:
        print("I failed to open a file")
    #if cr == None:
    #    raise IOError, "Unable to read in data from '%s'" % filename
    print >>f,"After Successful Read"
    ctype = get_crate_type(cr)
    if ctype != "Image":
        raise IOError, "The file '%s' does not contain an image (it is a %s)" % (filename,ctype)

    # Get the data.
    #
	#This needs to be get_piximgvals in CIAO 4.2 
    ivals = get_piximgvals(cr)
    dims = ivals.shape

    print 'd0 %d d1 %d d2 %d ' % (dims[0],dims[1],dims[2])

    #HOPS data is 3D so altering to 3D below
    # I guess I could let a 3D cube through...
    if len(dims) != 3:
        raise IOError, "Expected a 3D cube but found %d dims in '%s'" % (len(dims), filename)
    #if dims[0] != 1:
    #    raise IOError, "Expected the fourth dimension to have a length of 1, not %d ('%s')" % (dims[3],filename)
    
    nx = dims[1]
    ny = dims[2]
    nz = dims[0]

    # The reshape and copy are both to work around
    # Crates bugs in CIAO 4.0
    #
    
    #sl
    #ivals = ivals.reshape(nz,ny,nx).copy()
    #sl2
    ivals = ivals.reshape(nz,nx,ny).copy()
    
    #old
    #ivals = ivals.reshape(nx,ny,nz).copy()

    (v0,i0,dv) = get_velinfo(filename)
    vel = v0 + dv * ((numpy.arange(nz) + 1) - i0)

    return (ivals,vel)

def get_velinfo(filename):
    """Return the CRVAL, CRPIX, and CDELT values for the
    velocity axis of filename. To make things easy we assume
    that the velocity axis is the third axis.

    This is not very elegant
    """

    output = subprocess.Popen("dmlist '%s' opt=header,clean,raw" % filename,
                              shell=True, stdout=subprocess.PIPE).stdout
    cards = output.read().splitlines()

    # Inelegant Python warning
    #
    hdr = {}
    for c in cards:
        v = c.split("=",1)
        if len(v) == 2:
            hdr[v[0].strip()] = v[1].strip()

    def getval(str):
        "Extract the value of the keyword as a string - assumes no spaces"
        return str.split(" ",1)[0]


    def getfpval(str):
        "Extract the numeric (floating-point) value of the keyword - assumes no spaces"
        return float( getval(str) )


    if "CTYPE3" not in hdr:
        raise IOError, "Expected '%s' to contain a CTYPE3 keyword" % filename
    vname = getval(hdr["CTYPE3"])
    if vname != "VELO-LSR":
        raise IOError, "Expected CTYPE3 value to be VELO-LSR but found '%s'" % vname

    # I've gone functional-programming mad...
    return tuple (map( getfpval, (map (hdr.get, ["CRVAL3", "CRPIX3", "CDELT3"]))))

def write_image(filename,data,blockname="IMAGE"):
    """Create a FITS image file containing the data (>= 1d array).
    The blockname argument is used for the name of the block in the file
    (it defaults to IMAGE).

    This will clobber (i.e. overwrite) any existing file.
    """

    # We have to mangle the data so that the output is written correctly
    #
    cr = IMAGECrate()
    cd = CrateData()
    cd.name = blockname
    cd.set_nsets (1)

    dims = list (data.shape)
    print dims
    #dims.reverse()
    print dims
    f = open('junk.txt', 'a')

    print >>f,"I am a fish"
    print >>f,dims
    cd.set_dimarr(dims, len(dims))

    if add_piximg (cr, cd) != 1:
        raise ValueError, "Unable to add image block '%s' to a crate" % blockname

    if set_piximgvals(cr, data) != 1:
        raise ValueError, "Unable to add image (dimensions=%s) to a crate" % dims

    if write_file (cr, filename) != 1:
        raise IOError, "Unable to write image to file '%s'" % filename
