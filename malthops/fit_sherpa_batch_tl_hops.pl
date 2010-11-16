#! /usr/bin/perl -w
#$Id $
#$Log$
#------------------------------------------------------------
#usage:
#   ./fit_sherpa_batch_tl_hops.pl <INFILE>
#------------------------------------------------------------
#This program calls sherpa to do the fitting. It first creates
#an input file with all the fitting parameters in the correct
#format.
#------------------------------------------------------------

use strict;

 MAIN:
{  
    #------------------------------------------------------------
    #LOCAL VARS
    
    my ($name,$fits_file,$ngauss,$linkfwhm,$linkpos);
    my ($offsets,$maskimg,$plotname);
    
    #------------------------------------------------------------
    #READING INPUT FILE FROM CMD LINE AND SORTING FITTING PARS
    
    #Get command line args
    my $nfiles = scalar(@ARGV);
    if ($nfiles < 1) {
	die "\n\tNo input file given. Usage nh3_fit_sherpa_batch_tl.pl <INFILE> \n\n";
    }
    
    #Check file exists
    if (!(-e $ARGV[0])) {
	print "\'$ARGV[0] does not exist!\'\n";
	exit 0;
    }
    
    #Read input file
    open(IN, "<$ARGV[0]") || die "Can't open IN: $ARGV[0]\n";
    my @input =<IN>;
    
    #Searching for vals
    
    foreach(@input){
	unless(/^#/){
	    if(/^(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+$/){
		$name=$1;   $fits_file=$2;
		$ngauss=$3; $linkfwhm=$4; $linkpos=$5;
		$offsets=$6; $maskimg=$7; $plotname=$8;
		print "$1 $2 $3 $4 $5 $6 $7 $8\n";
	    }
	    else{
		die "\n\tERROR: Format not correct. Exiting.\n\n";
	    }
	

	    #------------------------------
	    #Checking format of each parameter is ok (do later)
	    
	    #------------------------------------------------------------
	    #RUNNING FITTING
	    
	    #creating sherpa scripts
	    open(SHERPA_FILE, ">sherpa_fit.py") || die "Can't open file.";
	    print SHERPA_FILE "#!/usr/bin/env python \n";

	    #snl no need for these loops as there will be no mask files (remove later)   
	    #if mask is fits file
	    #if( ($maskimg ne "None") && (!($maskimg =~ /\d/)) ){
	    if(1){
		#snl no mask files  
		#snl if (!(-e $maskimg)) {
		#snl    die "Error with mask file: $maskimg does not exist. Exit.\n"
		#snl}
		
		print SHERPA_FILE "import sys \n";
		print SHERPA_FILE "sys.path.append(\'/usr/bin\') \n";
		print SHERPA_FILE "sys.path.append(\'/usr/lib/python25.zip\') \n";
		print SHERPA_FILE "sys.path.append(\'/usr/lib/python2.5\') \n";
		print SHERPA_FILE "sys.path.append(\'/usr/lib/python2.5/plat-linux2\') \n";
		print SHERPA_FILE "sys.path.append(\'/usr/lib/python2.5/lib-tk\') \n";
		print SHERPA_FILE "sys.path.append(\'/usr/lib/python2.5/lib-dynload\') \n";
		print SHERPA_FILE "sys.path.append(\'/usr/local/lib/python2.5/site-packages\') \n";
		print SHERPA_FILE "sys.path.append(\'/usr/lib/python2.5/site-packages\') \n";
		print SHERPA_FILE "sys.path.append(\'/usr/lib/python2.5/site-packages/Numeric\') \n";
		print SHERPA_FILE "sys.path.append(\'/usr/lib/python2.5/site-packages/gst-0.10\') \n";
		print SHERPA_FILE "sys.path.append(\'/var/lib/python-support/python2.5\') \n";
		print SHERPA_FILE "sys.path.append(\'/usr/lib/python2.5/site-packages/gtk-2.0\') \n";
		print SHERPA_FILE "sys.path.append(\'/var/lib/python-support/python2.5/gtk-2.0\') \n";
		print SHERPA_FILE "sys.path.append(\'/var/lib/python-support/python2.5/IPython/Extensions\') \n";
		#### jfoster -- unnecessary extra path?
		#print SHERPA_FILE "sys.path.append(\'/home/slongmor/.ipython\') \n";
		####
		print SHERPA_FILE "import pyfits as pf\n";
		print SHERPA_FILE "import numpy as np\n";
		print SHERPA_FILE "import qfit_tl_hops as qfit\n";
		
		#snl no mask file so no need to create mask array
		#snl print SHERPA_FILE "maskarray=pf.getdata(\"$maskimg\") \n";
		#snl print SHERPA_FILE "maskarray=np.array(maskarray)\n";
		#print SHERPA_FILE "maskarray= maskarray.transpose() \n";
		#print SHERPA_FILE "execfile(\"qfit_tl_hops.py\") \n";

		#snl print SHERPA_FILE "res = qfit.fitfile(\"$fits_file\", ngauss=$ngauss, verbose=True,linkfwhm=$linkfwhm,linkpos=$linkpos,offsets=$offsets,maskimg=maskarray,plotname=$plotname) \n";
		print SHERPA_FILE "res = qfit.fitfile(\"$fits_file\", ngauss=$ngauss, verbose=True,linkfwhm=$linkfwhm,linkpos=$linkpos,offsets=$offsets,maskimg=$maskimg,plotname=$plotname) \n";
	    }
	    #snl should never get to the following loop but commenting out to be sure   
	    else{
		#snl print SHERPA_FILE "import qfit_tl_hops as qfit\n";
		#print SHERPA_FILE "execfile(\"qfit_tl_hops.py\") \n";
		#snl print SHERPA_FILE "res = qfit.fitfile(\"$fits_file\", ngauss=$ngauss, verbose=True,linkfwhm=$linkfwhm,linkpos=$linkpos,offsets=$offsets,maskimg=$maskimg,plotname=$plotname) \n";
	    }

	    print SHERPA_FILE "qfit.savefit(res,\"$name\") \n";
	    close(SHERPA_FILE);
	    
	    #running sherpa scripts
	    `python sherpa_fit.py`;
	    
	    #removing/cleaning up sherpa scripts
	    #`rm sherpa_fit.py`;
	}
    }
}

