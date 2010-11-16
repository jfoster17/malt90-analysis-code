
# Example of using the Python version of Sherpa, CIAO 4.0 to fit
# multiple gaussians
# to x,y data.
#
# We do not include a background component, but this is not too
# hard. We also do not exclude any data ranges, which could be
# done either before sending the data to setup_data (ie just
# don't send in any regions we do not want to fit) or
# by using ignore/notice once the data has been loaded
# into Sherpa.
#

from sherpa.astro.ui import *
from pychips.hlui import *
from pycrates import *
#import pycrates

from read_cube_tl_hops import *
#sl
from pychips import *

import logging
import types
import numpy

#execfile ("read_cube_tl_hops.py")
#execfile ("read_cube.py")

def includepos(x,y,minval=1.0):
    """Returns True if a fit should be tried for this spectrum,
    False if it should be excluded from the fit.

    This is a simple check that returns
        sum(y) >= minval

    Can we do this in a more pythony/numpy way?
    """

    return y.sum() > minval

def fitfile(cubename,ngauss=5,verbose=False,
          linkfwhm=False,linkpos=False,
          offsets=[-19600,-7450,0,7450,19600],
          maskimg=None,plotname=None
          ):
    """Fit ngauss gaussians to the data cube cubename.

    The meaning of the other parameters are the same
    as for the fitcube routine.
    """
    verbose = True
    if verbose:
        print "Reading in cube from %s" % cubename
    #f = open('junk.txt', 'a')
    #print >>f,cubename

    (cube, vaxis) = read_cube(cubename)
    #print >>f,"Cube check"
    

    res = fitcube(cube, vaxis,
                  ngauss, verbose,
                  linkfwhm, linkpos,
                  offsets, maskimg, plotname)
    return res

def savefit(results,prefix):
    """Save the results - the return value from fitfile() or
    fitcube() as FITS images. The files created are
      prefix + '.ampl.fits'
      prefix + '.fwhm.fits'
      prefix + '.pos.fits'
      prefix + '.resid.fits'
    
    At present there is no coordinate information included
    in the files (which is not ideal)."""

    if len(results) != 4:
        raise ValueError, "Expected a tuple of 4 elements for the results argument"

    def saveit(name, res):
        ofile = prefix + "." + name + ".fits"
        write_image (ofile, res, blockname=name)
        print "Created: %s" % ofile

    map (saveit, ["ampl", "fwhm", "pos", "resid"], results)

def fitcube(cube, vaxis,
            ngauss=5,verbose=False,
            linkfwhm=False,linkpos=False,
            offsets=[-19600,-7450,0,7450,19600],
            maskimg=None, plotname=None
            ):
    """Fit ngauss gaussians to the data in cube, given
    in (velocity,y,x) order, using as the velocity
    axis the vaxis array.

    The verbose flag controls the amount of screen output.

    The offsets array should be a sequence whose length equals
    ngauss. It determines the initial offsets of the peaks from
    the highest peak. If linkpos is True then the peaks are
    forced to use this relative spacing, otherwise the peak positions
    are fit independently.

    If linkfwhm is True then the fwhm of all the gaussians is
    fixed to the same value.

    If maskimg is set to None then a fit is attempted at all spatial
    points. If set, then it should be one of:
      a) an image, with the same spatial dimensions as cube, whose
         values are 1 for the spectrum at that point to be fitted,
         0 if excluded.
      b) a scalar value, which gives the minimum signal level at which
         a fit is allowed (integrated over all velocity components)
      c) a routine which takes two arguments, the velocity coordinates
         and the spectral values, and returns True if the point is
         to be fitted, False otherwise.
    
    If plotname is not None then postscript plots of fit+residuals
    are made at each successful point and called
        plotname + ".<i>.<j>.ps"
    where i and j are the x and y coordinates in pixels, starting
    at 0

    The return value is a tuple containing:
      best fit amplitudes
      best fit FWHM
      best fit positions
      residuals about the best fit (data - model)

    The ampliture array has dimensions ngauss by ny by nx.
    The FWHM array has dimensions      ny by nx
      IF linkfwhm is true, otherwise  ngauss by ny by nx.
    The pos array has dimensions       ny by nx
      IF linkpos is true, otherwise   ngauss by ny by nx.
    The residual cube has the same dimensions as cube, namely
                                      ngauss by ny by nx
    
    If linkpos is true then the reported position refers to
    the FIRST gaussian that is fitted (ie its offset from the
    largest peak is given by the offsets[0] value). This
    could be changed.

    """

    if len(cube.shape) != 3:
        raise ValueError, "cube argument must be a 3D array"
    if len(vaxis.shape) != 1:
        raise ValueError, "vaxis argument must be a 1D array"
    if ngauss < 1:
        raise ValueError, "ngauss must be >= 1"
    if len(offsets) != ngauss:
        raise ValueError, "ngauss must match the length of the offsets array"

    (nz,ny,nx) = cube.shape
    if len(vaxis) != nz:
        raise ValueError, "vaxis does not match the first dimension of cube"

    # slightly complicated check of maskimg argument
    # The logic disallows a scalar value being given as a numpy.ndarray value
    #
    if maskimg != None:
        if callable(maskimg):
            try:
                testval = maskimg(vaxis, cube[:,0,0])
            except Exception:
                raise ValueError, "Unable to call routine sent in via maskimg argument"
            if not(isinstance(testval,(types.BooleanType, numpy.bool_, numpy.bool))):
                raise ValueError, "Return value of routine sent in via maskimg argument not a boolean"

        elif isinstance(maskimg,numpy.ndarray):
            mdims = maskimg.shape
            if len(mdims) != 2:
                raise ValueError, "maskimg parameter not sent a 2D array"
            if mdims != (ny,nx):
                raise ValueError, "maskimg array size incorrect, expected (%d,%d)" % (ny,nx)

        elif isinstance(maskimg,(types.FloatType,types.IntType)):
            tmpval = maskimg
            maskimg = lambda x,y: includepos(x,y,minval=tmpval)

            try:
                testval = maskimg(vaxis, cube[:,0,0])
            except Exception:
                raise ValueError, "Unable to call includepos with minval set to '%s'" % tmpval

        else:
            raise ValueError, "maskimg parameter must be None, a scalar, a callable object, or a NumPy array"

    if verbose:
        print "Fitting gaussians to data cube" 
        print "  - dimensions are %d x %d x %d" % (nx,ny,nz)
        print "  - fitting with %d gaussians." % ngauss
        print "  - linkfwhm = %s" % linkfwhm
        print "  - linkpos  = %s" % linkpos
        print "  - linkfwhm = %s" % offsets
        print "  - maskimg = %s" % type(maskimg)

    # see if we can swap the cube order to make it more
    # efficient: change from vel,y,x to y,x,vel order
    # - perhaps should just change to x,y,vel order as
    #   this is one less swap
    # - also, should the order of the output cubes be
    #   changed to match (and then converted back at
    #   the end of the loop)?
    #
    # This did not seem to speed things up, but perhaps
    # that's because most elements weren't being fit in
    # the test - ie if given the whole cube the comparison
    # may be different.
    ###cube = cube.swapaxes(0,2).swapaxes(0,1).copy()

    # Output images
    #    cube for amplitudes (ngauss,ny,nx)
    #    if linkfwhm
    #      image for fwhm (ny,nx)
    #    else
    #      cube for fwhm (ngauss,ny,nx)
    #    if linkpos
    #      image for pos (ny,nx) - position of first peak
    #    else
    #      cube for pos (ngauss,ny,nx)
    #    cube for residuals (nz,ny,nx)
    #
    # I am sure this can be a lot neater with some not-particularly-clever
    # numpy indexing code
    #
    cdims = (ngauss,ny,nx)
    idims = (ny,nx)
    ampl = numpy.zeros(cdims, cube.dtype)
    residuals = numpy.zeros(cube.shape, cube.dtype)
    if linkfwhm:
        fwhm = numpy.zeros(idims, cube.dtype)
    else:
        fwhm = numpy.zeros(cdims, cube.dtype)
    if linkpos:
        pos = numpy.zeros(idims, cube.dtype)
    else:
        pos = numpy.zeros(cdims, cube.dtype)
    ampl *= numpy.nan
    fwhm *= numpy.nan
    pos *= numpy.nan
    residuals *= numpy.nan
    if verbose:
        print "  - created storage space"

    # Do we need to ensure plotting is done in "batch" mode?
    # We should restore this after the fit but do not bother
    # for now
    #
    if plotname != None:
        set_preference("window.display", "false")

    # Ugly code.
    # Is this the fastest way to loop through the data?
    #
    for y in range(ny):
        for x in range(nx):

            spec = cube[:,y,x]
            #
            # for use when swap the axis order of the cube
            #spec = cube[y,x,:]

            if callable(maskimg):
                mval = maskimg(vaxis,spec)
            elif isinstance(maskimg, numpy.ndarray):
                mval = maskimg[y,x] == 1
            else:
                mval = True

            if verbose:
                print "  - x/y = %d %d  mask = %s" % (x,y,mval) 

            if not(mval):
                continue

            # ugly bit of code to check the fit succeeded,
	    #
            flag = False
            try:
                res = fit_ngauss(vaxis, spec,
                                 num=ngauss, verbose=verbose,
                                 offsets=offsets,
                                 linkfwhm=linkfwhm, linkpos=linkpos)
                flag = True

            except Exception:
                # essentially a no-op
                # (we assume we can essentially ignore any
                # exception here, which isn't ideal)
                spec = None


            if not(flag):
		continue

	    if (not(res[3].succeeded)):
		continue

            # casting is probably a bit excessive here, but I want to
            # be explicit about what is going on
            #
            ampl[:,y,x] = numpy.asarray (res[2], dtype=cube.dtype)

            if linkfwhm:
                fwhm[y,x] = numpy.asarray (res[1][0], dtype=cube.dtype)
            else:
                fwhm[:,y,x] = numpy.asarray (res[1], dtype=cube.dtype)

            if linkpos:
                pos[y,x] = numpy.asarray (res[0][0], dtype=cube.dtype)
            else:
                pos[:,y,x] = numpy.asarray (res[0], dtype=cube.dtype)

            residuals[:,y,x] = numpy.asarray (res[4], dtype=cube.dtype)

	    if plotname == None:
	        continue

	    plot_fit_resid()
	    current_plot("plot1")
	    set_plot_title ("xpix = %d  ypix = %d" % (x,y))
	    set_plot_ylabel("Jy/beam")
	    current_plot("plot2")
	    set_plot_xlabel("Velocity (m s^{-1})")
	    set_plot_ylabel("Residuals")
	    yl = min(res[4])
	    yh = max(res[4])
	    yd = 0.05 * (yh - yl)
            limits (Y_AXIS, yl-yd, yh+yd)
            set_curve(["err.down",False,"err.up",False,
	    	      "symbol.style","none","line.style","solid"
	    	      ])
	    print_window("%s.%d.%d" % (plotname,x,y),
	    	      	["format","ps","fittopage",True,
	    	     	"orientation","landscape",
	    	     	"keepaspect",True,"pagesize","letter",
			"clobber",True])	


    return (ampl, fwhm, pos, residuals)

def find_peak(x,y):
    """Very simple estimator of the position, width, and normalisation
    of the largest peak in the data.

    At present it is not designed to work on noisy data, for
    datasets where the most-massive peak is 'confused', or where
    the peak is close to the edge of the data.

    Returns a tuple containing
      (pos, fwhm, ampl)
    or throws a ValueError if there is a problem.
    """

    # Find the position of the peak and estimate its FWHM
    # by finding the first points which fall below half the
    # peak's height. This assumes that the background is flat
    # and 0.
    #
    ymax = y.max()
    if (y==ymax).sum() > 1:
        print "Warning: Found multiple peaks; using the first one"
    idx = y.argmax()

    xmax = x[idx]
    hmax = ymax / 2.0

    xidx = (y <= hmax) & (x > xmax)
    if xidx.any() == False:
        raise ValueError, "find_peak unable to find upper limit to peak at x,y= %g %g" % (xmax,ymax)
    xhi = x[xidx][0]

    xidx = (y <= hmax) & (x < xmax)
    if xidx.any() == False:
        raise ValueError, "find_peak unable to find lower limit to peak at x,y= %g %g" % (xmax,ymax)
    xlo = x[xidx][-1]

    fwhm = xhi - xlo

    # We estimate the amplitude of the normalized gaussian using
    # the peak width and maximum value
    #
    ampl = ymax * fwhm * numpy.sqrt(0.125 * numpy.pi / numpy.log(2.0))

    return (xmax, fwhm, ampl)

def setup_plot_styles():
    "Sets up Sherpa to use our desired plot styles"

    # Missing fit plot preferences
    #
    for p in [get_data_plot_prefs(), get_model_plot_prefs()]:
        p["yerrorbars"] = False
        p["symbolstyle"] = chips_none
        p["linestyle"] = chips_solid

def setup_data(x,y):
    """Given two 1d arrays that represent a spectrum (x and y),
    set up Sherpa for fitting the data. We use the first dataset
    to load the data and set up fake errors of 1.0 for each bin.

    We do NOT set up/clear any models."""

    load_arrays (1, x, y)
    d = get_data (1)
    d.staterror = y * 0.0 + 1.0

def setup_model(names):
    """Create a set of normalized gaussians, with
    component names equal to the members of the names
    input argument.

    We set the source model for dataset 1 to be the
    sum of all these components.

    We assume that dataset 1 contains valid data as
    we set the min/max values for the position of
    the models based on the X axis range. The position
    of the components is evenly spread across this
    range.
    """

    d = get_data (1)
    x = d.x
    xlo = x.min()
    xhi = x.max()
    xsep = (xhi - xlo) / (1.0+len(names))
    xpos = xlo + xsep

    for n in names:
        create_model_component ("normgauss1d", n)
        p = get_par (n + ".pos")
        p.min = xlo
        p.max = xhi
        p.val = xpos
        xpos += xsep

    set_source (1, "+".join(names))


def adjust_ngauss_params(names,pos,offsets,fwhm,ampl,
                         linkfwhm=False,linkpos=False):
    """Given a set of normgauss1d model components whose
    names are given by the names array, set the
    position to pos+offsets (pos is a scalar and
    offsets a numpy array), the fwhm to fwhm
    (scalar) and the amplitute to ampl (scalar).

    If the optional linkfwhm argument is set to True
    (the default is False) then all the fwhm values
    are linked together.

    If the optional linkpos argument is set to True
    (the default is False) then the positions are
    fixed relative to each other, with the offsets
    defined by the offsets array, but there absolute
    positions can be shifted.

    We also adjust the min/max values for fwhm and
    ampl to */ 1e3.
    """

    # must be neater ways of doing this
    apos = pos + offsets
    afwhm = offsets * 0.0 + fwhm
    aampl = offsets * 0.0 + ampl

    # It does not matter which model is the one that
    # can vary if we link the parameters, so I
    # chose the first as it is the simplest solution.
    #
    # This was a nice simple loop but quickly got ugly
    #
    n0 = names[0]
    off0 = offsets[0]
    for (n,p,f,a,off) in zip(names,apos,afwhm,aampl,offsets):
        set_par(n + ".ampl", a, a/1.0e3, a*1.0e3)

        if linkpos and n != n0:
            link (n + ".pos", n0 + ".pos + %g" % (off-off0) )
        else:
            set_par(n + ".pos", p)

        if linkfwhm and n != n0:
            link (n + ".fwhm", n0 + ".fwhm")
        else:
            set_par(n + ".fwhm", f, f/1.0e3, f*1.0e3)

def get_ngauss_params(names):
    """Returns a tuple containing the position, fwhm, and amplitude
    values for the components with names given by the names array.

    There is no indication of whether the parameters were linked,
    or how they were linked.
    
    There is essentially no error checking.
    """

    res = [(get_par(n+".pos").val, get_par(n+".fwhm").val, get_par(n+".ampl").val)
           for n in names]
    # From
    #    http://paddy3118.blogspot.com/2007/02/unzip-un-needed-in-python.html?showComment=1181409120000#c3787283431220917975
    #
    return tuple(map(list,zip(*res)))

   
def fit_ngauss(x,y,num=5,verbose=False,
               linkfwhm=False,linkpos=False,
               offsets=[-19600,-7450,0,7450,19600]):
    """Fit num normalized gaussians to the (x,y) data
    values.

    If linkfwhm is set then the FWHM values are
    linked together (i.e. forced to be the same for
    each component).

    If linkpos is set then the relative positions of
    the gaussians are fixed at the values given by the
    offsets argument, although the absolute positions
    are allowed to vary.

    If linkpos is set to False then the offsets array is
    used to set up the initial positions of the gaussians
    relative to the brightest peak.

    If verbose is set then messages are printed out during
    the fit.

    The return value is a tuple containing
       a list of positions
       a list of FWHM values
       a list of amplitudes
       the fit results (the return value of Sherpa's get_fit_results)
       an array of residuals to the fit

    To see if the fit was successful (not that it necessarily produced
    a good fit, just that it finished), check the succeeded field of
    the fourth element of the tuple - it should be True.

    A ValueError will be raised if there was a problem in identifying the
    main peak (which occurs when setting up the initial parameter
    guesses).
    """

    if num < 1:
        raise ValueError, "num argument must be >= 1."
    if len(offsets) != num:
        raise ValueError, "num and offsets arguments do not match: num=%d len(offsets)=%d" % (num,len(offsets))

    # Not convinced this is sensible
    offsets = numpy.array(offsets)

    names = ["g%d" % (i+1) for i in range(num)]
    if verbose:
        print "Fitting %d normalized gaussians" % num

    (ipos, ifwhm, iampl) = find_peak(x,y)
    if verbose:
        print "Initial peak found: pos= %g  fwhm= %g  ampl= %g" % (ipos, ifwhm, iampl)

    setup_plot_styles()
    setup_data(x,y)
    if verbose:
        print "Loaded data"
    setup_model(names)
    if verbose:
        print "Normalized gaussian components: %s" % " ".join(names)
    adjust_ngauss_params(names,ipos,offsets,ifwhm,iampl,
                         linkfwhm=linkfwhm,linkpos=linkpos)
    if verbose:
        print "Adjusted initial parameters: linkfwhm=%s  linkpos=%s" % (linkfwhm,linkpos)

    if verbose:
        print "Starting fit"
        fit(1)
    else:
        quietfit(1)

    (pos,fwhm,ampl) = get_ngauss_params(names)
    return (pos, fwhm, ampl, get_fit_results(), get_ngauss_resid())

def quietfit(id=1):
    "Calls fit() but hides the screen output."

    slog = logging.getLogger("sherpa")
    oldval = slog.level
    slog.setLevel(0)
    fit(id)
    slog.setLevel(oldval)

def get_ngauss_resid():
    "Returns the residuals from the current fit to dataset 1."

    d = get_data(1)
    mdl = get_model(1)
    return d.y - mdl(d.x)

