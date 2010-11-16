#! /usr/bin/perl -w
#$Id $
#$Log$
#------------------------------------------------------------

#------------------------------------------------------------

#============================================================
#PREAMBLE

#reading in PDL and pgplot dirs
#### jfoster - Unnecessary?
BEGIN{push (@INC, "/mako4/b/py-extern/lib/perl/5.8.8/PGPLOT.pm");}
#BEGIN{push (@INC, "//PGPLOT-2.19/blib/lib/");}
BEGIN{push (@INC, "/usr/lib/perl5/PDL/Graphics/");}
BEGIN{push (@INC, "/usr/lib/perl5/PDL/Graphics/PGPLOT");}
BEGIN{push (@INC, "/usr/lib/perl5/PDL/");}
BEGIN{push (@INC, "/usr/lib/perl5/");}

#loading standard modules
use strict;
use Getopt::Long;
use PGPLOT;

#loading PDL modules
use PDL;
use PDL::Core;
use PDL::IO::FITS;
use PDL::AutoLoader;
use PDL::Slices;
use PDL::Graphics::PGPLOT;
use PDL::Graphics::PGPLOT::Window;
use PDL::Fit::Gaussian;
use PDL::Graphics::LUT;

#============================================================
 MAIN:
{
    open(LOG, ">log_file.txt") || die "Can't open file."; #opening log file
    
    #============================================================
    #DECLARING VARIABLES

    my ($nfiles,$latex_name, $transition);
    my ($source_name,$sherpa_name,$amp_K_file,$fwhm_kms_file,$vel_kms_file,$resid_K_file,$header_file);
    my ($hdr,$amp_K,$fwhm_kms,$vel_kms,$resid_K);
    my ($mean_ramp_fullcube,$prms_ramp_fullcube,$median_ramp_fullcube);
    my ($min_ramp_fullcube,$max_ramp_fullcube,$adev_ramp_fullcube,$rms_ramp_fullcube);
    my ($amp_K_sn_fullcube,@amp_K_dims,@fwhm_kms_dims, @resid_K_dims );
    my ($fwhm_tol_max_kms,$fwhm_tol_min_kms,$vlsr_tol_kms,$amp_tol_max_K,$amp_tol_min_K);
    my ($bad_fit_pix, $pix_change, $tmp, $bad_fit_mask);
    my ($amp_K_bfm,$fwhm_kms_bfm,$vel_kms_bfm,$resid_K_bfm);
    my ($at_least_one_good_fit,$map_or_not,$good_fit_pixs);
    my ($max_pix_amp_K_bfm, $nelem_good_fit_pixs,$min_pix_to_map);
    my ($bad_fit_mask_vlsr,$bad_fit_pix_vlsr,$good_analysis_at_max_amp);
    my ($mean_samp,$prms_samp,$median_samp,$min_samp,$max_samp,$adev_samp,$rms_samp);
    my ($mean_ramp,$prms_ramp,$median_ramp,$min_ramp,$max_ramp,$adev_ramp,$rms_ramp);

    my ($spectra_ps_out_dev,$spectra_ps_out_name);
    my ($spect_K_at_max_pos_bfm, $resid_K_at_max_pos_bfm);
    my ($max_pos_data, $spect_at_max_pos_data_K, $spect_to_plot, $resid_to_plot, $plot_resid);
    my ($spect_min_val, $res_min_val, $res_max_val, $resid_to_plot_dcoff, $spect_to_plot_vel_axis_kms);
    my (@spect_to_plot_vel_axis_kms_array, @spect_to_plot_array, @resid_to_plot_dcoff_array);
    my (@vel_kms_at_max_amp_bfm_array);
    my ($sum_resid, $sum_spect, @sep_line_x, @sep_line_y, @max_min);

    my ($col_lookup, $amp_name_out_bfm, $vel_name_out_bfm, $fwhm_name_out_bfm, $resid_name_out_bfm);
    my ($mean_amp,$prms_amp,$median_amp,$min_amp,$max_amp,$adev_amp,$rms_amp);
    my ($mean_vel,$prms_vel,$median_vel,$min_vel,$max_vel,$adev_vel,$rms_vel);
    my ($mean_fwhm,$prms_fwhm,$median_fwhm,$min_fwhm,$max_fwhm,$adev_fwhm,$rms_fwhm);
    my ($amp_miriad_name,$vel_miriad_name,$fwhm_miriad_name);
    my ($amp_K_at_max_amp_bfm,$vel_kms_at_max_amp_bfm,$fwhm_kms_at_max_amp_bfm);
    my ($output_name,$label,$caption);

    #============================================================
    #READING INPUT FILE FROM CMD LINE AND SORTING FITTING PARS
 
    #------------------------------
    #Get command line args
    $nfiles = scalar(@ARGV);
    if ($nfiles !=4) {
	die "\n\tIncorrect Input. Usage fit_sherpa_out_malt.pl <SOURCENAME> <SHERPA_NAME> <HEADER_FILE> <TRANSITION> \n\n"; 
    }
    
    $source_name   = $ARGV[0];
    $sherpa_name   = $ARGV[1];
    $amp_K_file    = $ARGV[1] . ".ampl.fits";
    $fwhm_kms_file = $ARGV[1] . ".fwhm.fits";
    $vel_kms_file  = $ARGV[1] . ".pos.fits";
    $resid_K_file  = $ARGV[1] . ".resid.fits";
    $header_file   = $ARGV[2];
    $transition    = $ARGV[3];

    #------------------------------
    #Check files exist
    if ( (!(-e $amp_K_file    )) || 
	 (!(-e $fwhm_kms_file )) || 
	 (!(-e $vel_kms_file  )) || 
	 (!(-e $resid_K_file  )) || 
	 (!(-e $header_file   )) 
	) {
	die "Error reading files. Exiting. \n";
    }

    #------------------------------
    #reading in and checking files

    $hdr      = rfits($header_file,"hdrcopy");
    $amp_K    = rfits($amp_K_file,"hdrcopy");     
    $fwhm_kms = rfits($fwhm_kms_file,"hdrcopy");  
    $vel_kms  = rfits($vel_kms_file,"hdrcopy");                            
    $resid_K  = rfits($resid_K_file,"hdrcopy");        
    
    #cutting down to 2 dims
    $amp_K     = $amp_K->slice(":,:,(0)");
    $fwhm_kms  = $fwhm_kms->slice(":,:,(0)");
    $vel_kms   = $vel_kms->slice(":,:,(0)");

    #converting to correct units
    $amp_K    *= 0.001;          #converting from [mK] to [K]
    $fwhm_kms *= 0.001;          #convert from [m/s] to [[km/s]
    $vel_kms  *= 0.001;          #convert from [m/s] to [[km/s]

    #------------------------------
    #getting rms -> using full residual cubes
    ($mean_ramp_fullcube,$prms_ramp_fullcube,$median_ramp_fullcube,
	$min_ramp_fullcube,$max_ramp_fullcube,$adev_ramp_fullcube,$rms_ramp_fullcube) 
	= stats($resid_K);

    print LOG "fullcube $mean_ramp_fullcube,$prms_ramp_fullcube,$median_ramp_fullcube,
               $min_ramp_fullcube,$max_ramp_fullcube,$adev_ramp_fullcube,$rms_ramp_fullcube \n";

    #outputting blind s/n map
    $amp_K_sn_fullcube = $amp_K/$rms_ramp_fullcube;
    $amp_K_sn_fullcube->wfits("amp_K_sn_fullcube.fits",-32);

    #checking dimensions of cubes
    @fwhm_kms_dims   = dims($fwhm_kms);    #
    @amp_K_dims      = dims($amp_K);       #
    @resid_K_dims    = dims($resid_K);

    print "fwhm_kms_dims = @fwhm_kms_dims\n";
    print "amp_K_dims    = @amp_K_dims \n";
    print "resid_K_dims  = @resid_K_dims \n";

    #------------------------------
    #setting tolerances    
    $fwhm_tol_max_kms = 10.0;                      #no fits greater than 10kms
    $fwhm_tol_min_kms = $hdr->hdr->{CRPIX3}/1000;  #channel width 
    $vlsr_tol_kms     = 5.0;                       #Vlsr tol - after get peak spect
    $amp_tol_max_K    = 5*max($hdr);               #dont accept fits >5*max in cube
    $amp_tol_min_K    = 5*$rms_ramp_fullcube;      #min is 5*rms
    $min_pix_to_map   = 30;                        #need at least double the beam area 

    print LOG "fwhm_tol_max_kms = $fwhm_tol_max_kms, fwhm_tol_min_kms = $fwhm_tol_min_kms \n";
    print LOG "amp_tol_max_K     = $amp_tol_max_K,     amp_tol_min_K    = $amp_tol_min_K \n";

    #------------------------------
    #finding bad pixels based on the tolerances

    $bad_fit_pix = whichND( ($fwhm_kms < $fwhm_tol_min_kms) | 
			    ($fwhm_kms > $fwhm_tol_max_kms) | 
			    ($amp_K < $amp_tol_min_K) | 
			    ($amp_K > $amp_tol_max_K) #| 
			    #($vel_kms < ($vel_at_max_amp_kms_bfm - $vlsr_tol_kms) ) |
			    #($vel_kms > ($vel_at_max_amp_kms_bfm + $vlsr_tol_kms) )
			    );

    print LOG "bad_fit_pix = $bad_fit_pix \n";

    #------------------------------
    #creating bad pixel masks

    #initialising to correct array size with all values equal to one
    $bad_fit_mask = ($amp_K * 0) +1;
    $bad_fit_mask->wfits("bad_fit_mask_init.fits",-32);

    #now changing to bad pixels 
    $pix_change=0;
    unless($bad_fit_pix->isempty){	
	$tmp = $bad_fit_mask->indexND($bad_fit_pix);
	$tmp .=0;
	#$bad_fit_mask = $bad_fit_mask->setbadif($bad_fit_mask==0);
	$bad_fit_mask->wfits("bad_fit_mask_update.fits",-32);
	#print LOG "changed pix mask \n $bad_fit_mask \n";
	$pix_change++;
    }

    if($pix_change==0) {
	$bad_fit_mask->wfits("bad_fit_msk_noupdate.fits",-32);
	print LOG "no pixels changed  \n";
    }

    #============================================================
    #APPLYING BAD FIT MASKS TO VARIABLES

    $amp_K_bfm        = $amp_K      * $bad_fit_mask;
    $fwhm_kms_bfm     = $fwhm_kms   * $bad_fit_mask;
    $vel_kms_bfm      = $vel_kms    * $bad_fit_mask;
    $resid_K_bfm      = $resid_K    * $bad_fit_mask;

    #============================================================
    #DETERMINE WHICH ANALYSIS/IMAGING REGIME WE ARE IN FOR THIS SOURCE
    #BASED ON THE BAD FIT MASK PARAMETERS CALCULATED ABOVE

    #checking for at least one good fit?
    $at_least_one_good_fit=$map_or_not=0;        #initialise with no good fits and no mapping
    $good_fit_pixs = whichND($bad_fit_mask==1);  #find any good pixs
    unless($good_fit_pixs->isempty){$at_least_one_good_fit = 1;} 

    #running loop for sources where there is at least one good fit
    if($at_least_one_good_fit==1){    

	#get pixel with maximum realiable fit amplitude
	$max_pix_amp_K_bfm = whichND( $amp_K_bfm==max($amp_K_bfm) );
	print LOG "max_pix_amp_K_bfm=$max_pix_amp_K_bfm \n";

	#enough pixels to map emission? check if num good pixs sig larger than beam area
	$nelem_good_fit_pixs = nelem($good_fit_pixs); #number of good fit elements (*2)
	if( $nelem_good_fit_pixs > $min_pix_to_map){   #enough good pixels to do mapping
	    $map_or_not = 1;
	    print LOG "map_or_not=$map_or_not; $nelem_good_fit_pixs > $min_pix_to_map \n"
	}
	else{
	    $map_or_not = 0;
	    print LOG "map_or_not=$map_or_not; $nelem_good_fit_pixs/2 < $min_pix_to_map \n"}

    }

    print LOG "at_least_one_good_fit = $at_least_one_good_fit \n";

    #============================================================
    #NOW PROCEED DEPENDING ON WHICH PARTICULAR REGIME THIS SOURCE IS IN

    $good_analysis_at_max_amp = 0.0;

    #------------------------------------------------------------
    #REGIME 1: NO DETECTION
    if($at_least_one_good_fit==0){ print LOG "No detection.\n"; }

    #------------------------------------------------------------
    #REGIME 2: DETECTION
    elsif($at_least_one_good_fit==1){ #DETECTION

	print LOG "Detection in at least one pixel. \n";

	#------------------------------
	#finding location of peak emission and get spectra + fit resids 
	$max_pix_amp_K_bfm       = whichND($amp_K_bfm==max($amp_K_bfm));
	$spect_K_at_max_pos_bfm  = $hdr($max_pix_amp_K_bfm->at(0,0), 
						    $max_pix_amp_K_bfm->at(1,0),:)->sever;   
	$resid_K_at_max_pos_bfm  = $resid_K($max_pix_amp_K_bfm->at(0,0),
							    $max_pix_amp_K_bfm->at(1,0),:)->sever;   

	$amp_K_at_max_amp_bfm    = $amp_K_bfm->at(  $max_pix_amp_K_bfm->at(0,0), 
						 		$max_pix_amp_K_bfm->at(1,0));
	$vel_kms_at_max_amp_bfm  = $vel_kms_bfm->at($max_pix_amp_K_bfm->at(0,0), 
								$max_pix_amp_K_bfm->at(1,0));
	$fwhm_kms_at_max_amp_bfm = $fwhm_kms_bfm->at($max_pix_amp_K_bfm->at(0,0), 
								 $max_pix_amp_K_bfm->at(1,0));

	print LOG "$max_pix_amp_K_bfm \n";
	#print LOG "spect_K_at_max_pos_bfm = $spect_K_at_max_pos_bfm \n";
	#print LOG "resid_K_at_max_pos_bfm = $resid_K_at_max_pos_bfm \n";
	print LOG "vel_kms_at_max_amp_bfm  = $vel_kms_at_max_amp_bfm\n";
	print LOG "fwhm_kms_at_max_amp_bfm = $fwhm_kms_at_max_amp_bfm \n";
	
	#------------------------------
	#updating bad fit mask for fits with large vlsr offset from peak vlsr
	$bad_fit_mask_vlsr = $bad_fit_mask; 
 
	$bad_fit_pix_vlsr = whichND( ($vel_kms_bfm   < ($vel_kms_at_max_amp_bfm   - $vlsr_tol_kms) ) |
					($vel_kms_bfm   > ($vel_kms_at_max_amp_bfm   + $vlsr_tol_kms) ) );
	
	
	print LOG "$vel_kms_at_max_amp_bfm   - $vlsr_tol_kms \n";
	print LOG "$vel_kms_at_max_amp_bfm   + $vlsr_tol_kms\n";

	unless($bad_fit_pix_vlsr->isempty){	
	    $tmp = $bad_fit_mask_vlsr->indexND($bad_fit_pix_vlsr);
	    $tmp .=0;
	    $bad_fit_mask_vlsr->wfits("bad_fit_msk_vlsr.fits",-32);
	    $pix_change++;
	}
	
	#------------------------------
	#applying vlsr cut
		
	$amp_K_bfm        = $amp_K_bfm      * $bad_fit_mask_vlsr;
	$amp_K_bfm        = $amp_K_bfm * $bad_fit_mask_vlsr;
	$fwhm_kms_bfm     = $fwhm_kms_bfm   * $bad_fit_mask_vlsr;
	$vel_kms_bfm      = $vel_kms_bfm    * $bad_fit_mask_vlsr;
	$resid_K_bfm      = $resid_K_bfm    * $bad_fit_mask_vlsr;
	
	#------------------------------
        #getting statistics of spectra and resids
	($mean_samp,$prms_samp,$median_samp,$min_samp,$max_samp,$adev_samp,$rms_samp) 
	    = stats($spect_K_at_max_pos_bfm);
	($mean_ramp,$prms_ramp,$median_ramp,$min_ramp,$max_ramp,$adev_ramp,$rms_ramp) 
	    = stats($resid_K_at_max_pos_bfm);
	$sum_resid=sum($resid_K_at_max_pos_bfm);
	$sum_spect=sum($spect_K_at_max_pos_bfm);
	print LOG "peak spect $mean_samp,$prms_samp,$median_samp,$min_samp,$max_samp,
               $adev_samp,$rms_samp \n";
	print LOG "peak resid $mean_ramp,$prms_ramp,$median_ramp,$min_ramp,$max_ramp,
               $adev_ramp,$rms_ramp \n";
	open(STATS, ">>stats.txt") || die "Can't open file."; #outputting to stats file
	print STATS "$rms_ramp \t $sum_resid \t $sum_spect \t ";
	print STATS "$sum_resid \t $sum_spect \t\n";
	close(STATS);

    }

    #============================================================
    #NOW HAVE DERIVED VALUES, UPPER LIMITS OR NON DETECTIONS FOR EACH
    #PARAMETER, OUTPUT TO THE LATEXTABLE FILE

    open (LATEXTABLE, ">>fits.tex") || die "Can't open latex file.";

    #getting rid of underscores in name as latex doesn't like this
    $latex_name = $source_name;
    $latex_name =~ s/_//g; 
    print LATEXTABLE "$latex_name ";

    #------------------------------
    #non-detection: dash for everything
    if($at_least_one_good_fit==0){ 
	
	#                    TA    Vlsr  DV  RMS
	print LATEXTABLE " & -   & -   & -   & -   "; 
    }

    #------------------------------
    #detection
    else{ 
	#                     TA      Vlsr    DV     RMS
	printf LATEXTABLE " & %4.2f & %4.2f & %4.2f & %4.2f ", $amp_K_at_max_amp_bfm,
                                                               $vel_kms_at_max_amp_bfm,
                                                               $fwhm_kms_at_max_amp_bfm,
	                                                       $rms_ramp_fullcube;
    }

    print LATEXTABLE " \\\\ \n";    
    close(LATEXTABLE);

    #============================================================
    #CREATING POSTSCRIPT FILES OF PEAK SPECTRA

    #creating ps device file name and output ps filename
    $spectra_ps_out_dev = $source_name . "_spect.ps/ps";
    $spectra_ps_out_name = $source_name . "_spect.ps";

    #------------------------------
    #non-detection: no residual to fit
    if($at_least_one_good_fit == 0){

	$max_pos_data            = whichND($hdr==max($hdr));
	$spect_at_max_pos_data_K = $hdr($max_pos_data->at(0,0), $max_pos_data->at(1,0),:)->sever; 

	$spect_to_plot = $spect_at_max_pos_data_K;  #setting both resid and spect to same as
	$resid_to_plot = $spect_at_max_pos_data_K;  #will not plot resid
	$plot_resid    = 0;

	#print LOG "spect_to_plot=$spect_to_plot \n";
	#print LOG "resid_to_plot=$resid_to_plot \n";
    }

    #------------------------------
    #detection: 
    else{
	$plot_resid    = 1;
	$spect_to_plot = $spect_K_at_max_pos_bfm;
	$resid_to_plot = $resid_K_at_max_pos_bfm;
	
	#print LOG "spect_to_plot=$spect_to_plot \n";
	#print LOG "resid_to_plot=$resid_to_plot \n";
    }

    #getting axes ranges 
    $spect_min_val       = min($spect_to_plot); 
    $res_min_val         = min($resid_to_plot);
    $res_max_val         = max($resid_to_plot);
    $resid_to_plot_dcoff = $resid_to_plot -
	$plot_resid * ( sqrt($spect_min_val*$spect_min_val) +  sqrt($res_max_val*$res_max_val)   );

    #------------------------------
    #Calculating velocity of each pixel from header info
    $spect_to_plot_vel_axis_kms = 0.001 * (
	                           ( (sequence($hdr->hdr->{NAXIS3})
	                             - $hdr->hdr->{CRPIX3}) *
	                             $hdr->hdr->{CDELT3} ) +
	                           $hdr->hdr->{CRVAL3}
	                           );

    #------------------------------
    #Creating ps file

    #converting pdl to array so can be read into pgplot
    @spect_to_plot_vel_axis_kms_array = list $spect_to_plot_vel_axis_kms;
    @spect_to_plot_array              = list $spect_to_plot;
    @resid_to_plot_dcoff_array        = list $resid_to_plot_dcoff;

    @max_min = ( $spect_min_val-1.1*($res_max_val - $res_min_val),max($hdr));

    #creating an array with vlsr of peak emission
    if($at_least_one_good_fit == 0){
	print LOG "max_pos_data = $max_pos_data  \n";
    }
    else{
	@vel_kms_at_max_amp_bfm_array = ($vel_kms_at_max_amp_bfm,$vel_kms_at_max_amp_bfm);    
	print LOG "max_pix_amp_K_bfm = $max_pix_amp_K_bfm \n";
	print LOG "vel_kms_at_max_amp_bfm_array=@vel_kms_at_max_amp_bfm_array \n";
	    
    }

    print LOG "plot_resid = $plot_resid \n";
    print LOG "spect_min_val = $spect_min_val: res_min_val = $res_min_val: res_max_val=$res_max_val\n";
    #print LOG "spect_to_plot = $spect_to_plot \n";
    #print LOG "resid_to_plot = $resid_to_plot \n";
    #print LOG "resid_to_plot_dcoff = $resid_to_plot_dcoff \n";
    #print LOG "spect_to_plot_vel_axis_kms = $spect_to_plot_vel_axis_kms \n";
    #print LOG "spect_to_plot_array = @spect_to_plot_array \n";
    #print LOG "resid_to_plot_dcoff_array = @resid_to_plot_dcoff_array \n";
    #print LOG "max_min=@max_min \n";

    #------------------------------
    #creating figure

	# jbf -- set to stand-alone PGPLOT
	$ENV{'PGPLOT_DIR'} = '/usr/lib/pgplot5/';
	
    pgopen($spectra_ps_out_dev);
    pgpap(10,1.0);
    pgsvp(0.13,0.99,0.01,0.99);
    pgwnad(0.0,1.0,0.0,0.7);
    pgslw(3.0);  pgscf(1); pgsch(1.3);
    pgswin($spect_to_plot_vel_axis_kms->at(0),$spect_to_plot_vel_axis_kms->at($hdr->hdr->{NAXIS3} - 1),
	   ( $spect_min_val-1.1*$plot_resid*($res_max_val - $res_min_val)), max($hdr));
    pgtbox("bcnts",0.0,0.0,"vbncts",0.,0.);
    pgline($hdr->hdr->{NAXIS3},\@spect_to_plot_vel_axis_kms_array,\@spect_to_plot_array);
    pgline($hdr->hdr->{NAXIS3},\@spect_to_plot_vel_axis_kms_array,\@resid_to_plot_dcoff_array);
    pgsls(2);
    if($at_least_one_good_fit == 1){
	pgline(2, \@vel_kms_at_max_amp_bfm_array, \@max_min);
    }
    pgsls(4);
    if($plot_resid){ #plotting line to separate spectra frum resid
	@sep_line_x = ($spect_to_plot_vel_axis_kms->at(0), $spect_to_plot_vel_axis_kms->at($hdr->hdr->{NAXIS3}-1));
	@sep_line_y = ($spect_min_val,$spect_min_val);
	pgline(2, \@sep_line_x, \@sep_line_y);
	print LOG "sep_line_x=@sep_line_x sep_line_y=@sep_line_y\n";
    }
    
    pgmtxt("L",3.0,0.5,0.5,"T\\dA\\u\\u*\\d [K]");
    pgmtxt("B",3.0,0.5,0.5,"V\\dLSR\\u [kms\\u-1\\d]");
    pgmtext("T",-1.4,0.05,0.0,"$source_name");
    pgmtext("T",-1.4,0.95,1.0,"$transition");
    pgend();  
	# jbf -- back to MIRIAD PGPLOT
	$ENV{'PGPLOT_DIR'} = '/usr/local/miriad/lib/linux';
	

    #------------------------------
    #outputting appropriate info to latex file
    open (LATEXSPECTRA, ">>spectra_fig.tex") || die "Can't open file.";
    print LATEXSPECTRA "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$spectra_ps_out_name} \\\\  \n";
    close(LATEXSPECTRA);

    #============================================================
    #DECIDING WHETHER OR NOT TO MAP

    $col_lookup = -8;   #setting miriad colour lookup table number

    if($map_or_not == 1){ #yes, create maps

	#------------------------------
	#output fits where bad fit mask applied
	$amp_name_out_bfm   = $source_name . "_amp_bfm.fits"; 
	$vel_name_out_bfm   = $source_name . "_vel_bfm.fits"; 
	$fwhm_name_out_bfm  = $source_name . "_fwhm_bfm.fits"; 
	$resid_name_out_bfm = $source_name . "_resid_wcs_bfm.fits"; 
	
	badmask($amp_K_bfm->inplace,0);
	badmask($vel_kms_bfm->inplace,0);
	badmask($fwhm_kms_bfm->inplace,0);
	badmask($resid_K_bfm->inplace,0);

	$amp_K_bfm    = $amp_K_bfm->setbadif($amp_K_bfm==0);
	$vel_kms_bfm  = $vel_kms_bfm->setbadif($vel_kms_bfm==0);
	$fwhm_kms_bfm = $fwhm_kms_bfm->setbadif($fwhm_kms_bfm==0);
	$resid_K_bfm  = $resid_K_bfm->setbadif($resid_K_bfm==0);

	#copying headers
	cp_hdr($amp_K_bfm,    $hdr,"K");
	cp_hdr($vel_kms_bfm,  $hdr,"km/s");    
	cp_hdr($fwhm_kms_bfm, $hdr,"km/s");
	cp_hdr($resid_K_bfm,  $hdr,"K");
	
	$amp_K_bfm->wfits($amp_name_out_bfm,-32);		
	$vel_kms_bfm->wfits($vel_name_out_bfm,-32);
	$fwhm_kms_bfm->wfits($fwhm_name_out_bfm,-32);
	$resid_K_bfm->wfits($resid_name_out_bfm,-32);
	
	#------------------------------
	#getting limits
	($mean_amp,$prms_amp,$median_amp,$min_amp,$max_amp,$adev_amp,$rms_amp) 
	    = stats($amp_K_bfm);
	($mean_vel,$prms_vel,$median_vel,$min_vel,$max_vel,$adev_vel,$rms_vel) 
	    = stats($vel_kms_bfm);
	($mean_fwhm,$prms_fwhm,$median_fwhm,$min_fwhm,$max_fwhm,$adev_fwhm,$rms_fwhm) 
	    = stats($fwhm_kms_bfm);
	
	print LOG "min_amp  = $min_amp   max_amp  = $max_amp \n";
	print LOG "min_vel  = $min_vel   max_vel  = $max_vel\n";
	print LOG "min_fwhm = $min_fwhm  max_fwhm = $max_fwhm\n";
	
	#------------------------------
	#miriad steps
	$amp_miriad_name=$amp_name_out_bfm; 
	$amp_miriad_name =~s/fits/xy/;
	`fits in=$amp_name_out_bfm out=$amp_miriad_name op=xyin`;
	print LOG "in=$amp_name_out_bfm out=$amp_miriad_name op=xyin \n";
	
	$vel_miriad_name=$vel_name_out_bfm; 
	$vel_miriad_name =~s/fits/xy/;
	`fits in=$vel_name_out_bfm out=$vel_miriad_name op=xyin`;
	print LOG "in=$vel_name_out_bfm out=$vel_miriad_name op=xyin \n";
	
	$fwhm_miriad_name=$fwhm_name_out_bfm; 
	$fwhm_miriad_name =~s/fits/xy/;
	`fits in=$fwhm_name_out_bfm out=$fwhm_miriad_name op=xyin`;
	print LOG "in=$fwhm_name_out_bfm out=$fwhm_miriad_name op=xyin \n";
	
	#amp on amp
	 `cgdisp in=$amp_miriad_name,$amp_miriad_name type=p,c slev=p,1 levs1=10,20,30,40,50,60,70,80,90 device=$source_name.amp_amp.ps/cps  beamtyp=b,l,2  options=wedge,blacklab labtyp=absdeg,absdeg csize=1.4,0,0,0 range=$min_amp,$max_amp,lin,$col_lookup cols1=1`;


	#vel on amp colour
	`cgdisp in=$vel_miriad_name,$amp_miriad_name type=p,c slev=p,1 levs1=10,20,30,40,50,60,70,80,90 device=$source_name.vel_amp.ps/cps  beamtyp=b,l,2  options=wedge,blacklab  labtyp=absdeg,absdeg csize=1.4,0,0,0 range=$min_vel,$max_vel,lin,$col_lookup cols1=1`;
	
	#fwhm on amp
	`cgdisp in=$fwhm_miriad_name,$amp_miriad_name type=p,c slev=p,1 levs1=10,20,30,40,50,60,70,80,90 device=$source_name.fwhm_amp.ps/cps  beamtyp=b,l,2  options=wedge,blacklab  labtyp=absdeg,absdeg csize=1.4,0,0,0 range=$min_fwhm,$max_fwhm,lin,$col_lookup cols1=1`;

	#------------------------------
	#create latex file
	$output_name = $source_name . "output_fig.tex";
	$label="fig:" . $latex_name . "_fig";
	$caption = "\\small{$transition emission towards
                   source $latex_name. [Top left] 
                   $transition spectra with
                   residuals from the fit displayed underneath.
                   [Top right] Integrated intensity 
                   emission overlayed with contours in 10\\% steps of the peak
                   value, from 90\\% down. The bottom left and right panels
                   the FWHM
                   and V\$_{\\rm LSR}\$ of the emission determined
                   by the fit.} ";
                   
	open(OUTPUT, ">$output_name") || die "Can't open file.";
	print OUTPUT "\\documentclass[12pt,psfig]{mn2e} \n";
	print OUTPUT "\\usepackage{graphicx} \n";
	print OUTPUT "\\begin{document} \n";
	print OUTPUT "\\begin{figure*} \n";
	print OUTPUT "\\begin{center} \n";
	print OUTPUT "\\begin{tabular}{cc} \n";
	print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$spectra_ps_out_name} &  \n";
	print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source_name.amp_amp.ps}  \\\\ \n";
	print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source_name.fwhm_amp.ps} &  \n";
	print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source_name.vel_amp.ps} \\\\ \n";
	print OUTPUT "\\end{tabular} \n";
	print OUTPUT "\\end{center} \n";
	print OUTPUT "\\caption{$caption}\n";
	print OUTPUT "\\label{$label} \n";
	print OUTPUT "\\end{figure*} \n";
	print OUTPUT "\\end{document} \n";
	close(OUTPUT);
    }
    else{  #no mapping so only show spectra

	#------------------------------
	#create latex file
	$output_name = $source_name . "output_fig.tex";
	$label="fig:" . $latex_name . "_fig";
	$caption = "\\small{Non-detection of $transition emission towards
                   source $latex_name. } ";
                   
	open(OUTPUT, ">$output_name") || die "Can't open file.";
	print OUTPUT "\\documentclass[12pt,psfig]{mn2e} \n";
	print OUTPUT "\\usepackage{graphicx} \n";
	print OUTPUT "\\begin{document} \n";
	print OUTPUT "\\begin{figure*} \n";
	print OUTPUT "\\begin{center} \n";
	print OUTPUT "\\begin{tabular}{cc} \n";
	print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$spectra_ps_out_name} &  \n";
	#print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source_name.amp_amp.ps}  \\\\ \n";
	#print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source_name.fwhm_amp.ps} &  \n";
	#print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source_name.vel_amp.ps} \\\\ \n";
	print OUTPUT "\\end{tabular} \n";
	print OUTPUT "\\end{center} \n";
	print OUTPUT "\\caption{$caption}\n";
	print OUTPUT "\\label{$label} \n";
	print OUTPUT "\\end{figure*} \n";
	print OUTPUT "\\end{document} \n";
	close(OUTPUT);

    }
    
    #------------------------------
    #running latex to create composite postscript
    `latex $output_name`; 
    $output_name=~s/tex/dvi/;   
    `dvips $output_name`;
    `rm *output_fig.aux *output_fig.dvi *output_fig.log *output_fig.tex`;

    #tidying up
    if($map_or_not == 1){`rm -r $amp_miriad_name $vel_miriad_name $fwhm_miriad_name`;}
    
    close(LOG);
}


#------------------------------------------------------------
#SUBROUTINES

#adding header info
sub cp_hdr {
    my $fits_1 = shift;
    my $fits_2 = shift;
    my $bunit  = shift;
    #my $naxis  = shift;
    
    print LOG "bunit=$bunit\n";
    
    #print LOG "fits\n$fits_1 \n";
    #my $test=$fits_2->hdr->{CTYPE3};
    #print LOG "ctype3=$test \n";

    $fits_1->hdr->{CTYPE1} = $fits_2->hdr->{CTYPE1};
    $fits_1->hdr->{CTYPE2} = $fits_2->hdr->{CTYPE2};
    $fits_1->hdr->{CTYPE3} = $fits_2->hdr->{CTYPE3};
    $fits_1->hdr->{CRPIX1} = $fits_2->hdr->{CRPIX1};
    $fits_1->hdr->{CRPIX2} = $fits_2->hdr->{CRPIX2};
    $fits_1->hdr->{CRPIX3} = $fits_2->hdr->{CRPIX3};
    $fits_1->hdr->{CDELT1} = $fits_2->hdr->{CDELT1};
    $fits_1->hdr->{CDELT2} = $fits_2->hdr->{CDELT2};
    $fits_1->hdr->{CDELT3} = $fits_2->hdr->{CDELT3};
    $fits_1->hdr->{CRVAL1} = $fits_2->hdr->{CRVAL1};
    $fits_1->hdr->{CRVAL2} = $fits_2->hdr->{CRVAL2};
    $fits_1->hdr->{CRVAL3} = $fits_2->hdr->{CRVAL3};    
    $fits_1->hdr->{BPA}    = $fits_2->hdr->{BPA};
    $fits_1->hdr->{BMIN}   = $fits_2->hdr->{BMIN};
    $fits_1->hdr->{BMAJ}   = $fits_2->hdr->{BMAJ};
    $fits_1->hdr->{BUNIT}  = $bunit;

    #$fits_1->hdr->{EPOCH}=$fits_2->hdr->{EPOCH};
    #$fits_1->hdr->{BTYPE}=$fits_2->hdr->{BTYPE};
}
