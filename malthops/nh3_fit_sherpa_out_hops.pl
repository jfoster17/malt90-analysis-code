#! /usr/bin/perl -w
#$Id $
#$Log$
#------------------------------------------------------------
#./nh3_fit_sherpa_out.pl name 11_sherpout_name 22_sherpout_name 11_orig_file 22_orig_file mask_file_11 mask_file_22 fit_region
#
#e.g. ./nh3_fit_sherpa_out.pl hops_reg1 source_1_11 source_1_22 G330_a_NH3_clumps/NH311/clump_1_img.fits G330_a_NH3_clumps/NH322/clump_1_img.fits G330_a_NH3_clumps/NH311/clump_1_msk.fits G330_a_NH3_clumps/NH311/clump_1_msk.fits none

#./nh3_fit_sherpa_out_hops.pl G336.159+0.214_-151673.7_img G336.159+0.214_-151673.7_img_11 G336.159+0.214_-151673.7_img_22 /home/slongmor/sherpa/ciao-4.1/hops/crp_test1/G330_all_clumps/NH3_11/tmp2/G336.159+0.214_-151673.7_img_3.fits /home/slongmor/sherpa/ciao-4.1/hops/crp_test1/G330_all_clumps/NH3_22/tmp2/G336.159+0.214_-151673.7_img_5.fits  /home/slongmor/sherpa/ciao-4.1/hops/crp_test1/G330_all_clumps/NH3_11/tmp2/G336.159+0.214_-151673.7_msk.fits /home/slongmor/sherpa/ciao-4.1/hops/crp_test1/G330_all_clumps/NH3_22/tmp2/G336.159+0.214_-151673.7_msk.fits none
#------------------------------------------------------------
#To do
#
#------------------------------------------------------------

#============================================================
#PREAMBLE

#reading in PDL and pgplot dirs
#### jfoster - Unnecessary?
BEGIN{push (@INC, "/mako4/b/py-extern/lib/perl/5.8.8/PGPLOT.pm");}
#BEGIN{push (@INC, "//PGPLOT-2.19/blib/lib/");}
#BEGIN{push (@INC, "/home/slongmor/programming/PDL/PDL-2.4.2/Graphics/PGPLOT/");}
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
#PHYSICAL CONSTANTS AND NH3 PARAMTERS NEEDED FOR DERIVATIONS
#DECLARED AS GLOBALS

our $c_ms=2.99792458E8;         #[m/s]
our $pi=3.1415926535;           #[-] 
our $h=6.626068E-34;            #[m^2 kg / s]
our $k=1.3806503E-23;           #[m^2 kg s-2 K-1]

#------------------------------
#NH3 parameters
#converting 1.476 debye to SI = 1.476 * 3.33564E-30 Coulomb metres
our $nh3_dipole_moment_Cm = 4.92340464E-30; 

#freqncy of transitions in Hz (s^-1)
our @freq_Hz = (23.6945E9, 23.7226E9, 23.8701E9, 24.1394E9, 24.5330E9, 25.0560E9);

#Energy of the levels above ground expressed as a temperature in K 
our @T_above_ground_K = (22.71, 63.91, 123.03, 200.02, 294.89, 407.60);
our @T_diff_K;

our ($i, $j);
for($i=0;$i<6;$i++){
    for($j=0;$j<6;$j++){
	$T_diff_K[$i][$j]=$T_above_ground_K[$j]-$T_above_ground_K[$i];
    }
}

#theoretical fraction of the transition intensity in the main line
our @frac_main_sat = (0.5003, 0.796, 0.894, 0.935, 0.956, 0.969);

#einstein A co- efficients for the transitions (transition rate (s^-1)
our @Einstein_A_s = (1.71E-7, 2.29E-7, 2.62E-7, 2.89E-7, 3.16E-7, 3.46E-7); #[s^-1]

#degeneracy of the transitions 
our @g=(6, 10, 28, 18, 22, 52);

#============================================================
 MAIN:
{
    open(LOG, ">log_file.txt") || die "Can't open file."; #opening log file
    
    #============================================================
    #DECLARING VARIABLES

    #source names and files holding data + fit results
    my ($nfiles, $source, $source_only, $source_tex, $l_deg, $b_deg, $region, $tmp);
    my ($nh3_11_amp_K_file, $nh3_11_fwhm_kms_file, $nh3_11_vel_kms_file, $nh3_11_resid_K_file);
    my ($nh3_22_amp_K_file, $nh3_22_fwhm_kms_file, $nh3_22_vel_kms_file, $nh3_22_resid_K_file);
    my ($header_file_11, $header_file_22, $duchamp_mask_file_11, $duchamp_mask_file_22 );

    #variables (pdl's) to store data and fit results
    my ($hdr_11, $hdr_22);                #11 and 22 data files -> also used to extract correct fits header info
    my ($nh3_11_amp_K, $nh3_11_amp_main_K, $nh3_11_fwhm_kms, $nh3_11_vel_kms, $nh3_11_resid_K, $nh3_11_amp_rms_K);
    my ($nh3_22_amp_K, $nh3_22_amp_main_K, $nh3_22_fwhm_kms, $nh3_22_vel_kms, $nh3_22_resid_K, $nh3_22_amp_rms_K);
    my ($nh3_11_duchamp_mask, $nh3_22_duchamp_mask);  #masks from emission-finding routine -> selects pix with emission 

    #statistics of post-fit residuals for the whole 11 and 22 cubes 
    # - used to get rms amplitude for each transition 
    # - this rms is one of the criteria used to distinguish between reliable vs non-reliable fits.
    my ($mean_ramp11_fullcube,$prms_ramp11_fullcube,$median_ramp11_fullcube,
	$min_ramp11_fullcube,$max_ramp11_fullcube,$adev_ramp11_fullcube,$rms_ramp11_fullcube);
    my ($mean_ramp22_fullcube,$prms_ramp22_fullcube,$median_ramp22_fullcube,
	$min_ramp22_fullcube,$max_ramp22_fullcube,$adev_ramp22_fullcube,$rms_ramp22_fullcube); 
    my ($nh3_11_amp_main_K_sn_fullcube, $nh3_22_amp_main_K_sn_fullcube);

    #tolerance values used to determine 
    # - whether a fit is reliable or not
    # - whether or not there are enough pixels to do mapping
    my ($fwhm_tol_max_11_kms, $fwhm_tol_min_11_kms, $vlsr_tol_11_kms, $amp_tol_max_11_K, $amp_tol_min_11_K);
    my ($fwhm_tol_max_22_kms, $fwhm_tol_min_22_kms, $vlsr_tol_22_kms, $amp_tol_max_22_K, $amp_tol_min_22_K);
    my ($min_pix_to_map);    #minimum number of pixels required before create maps of parameters

    #variables used in separating pixels with reliable vs unreliable fits
    my (@nh3_11_fwhm_kms_dims, @nh3_11_amp_main_K_dims,                   #intermediate variables to check the dimensions 
	@nh3_22_fwhm_kms_dims, @nh3_22_amp_main_K_dims);                  # of the input data are all the same
    my ($bad_fit_pix_11,        $bad_fit_pix_22);                         #list of pix with unreliable fits based on tols above
    my ($good_fit_pix_sherpa_11,  $good_fit_pix_sherpa_22);  
    my ($bad_fit_mask_sherpa_11, $bad_fit_mask_sherpa_22);                #list of pix with bad fits after sherpa fitting
    my ($bad_fit_pix_11_vlsr, $bad_fit_pix_22_vlsr);                      #list of pix with unreliable fits based on vlsr tols above
    my ($pix_change_11, $pix_change_22);                                  # = 1 or 0 depending if additional pix flagged bad or not
    my ($bad_fit_mask_11,      $bad_fit_mask_22,      $bad_fit_mask_11_22_comb);         #masks made from bad pixel lists 
    my ($bad_fit_mask_11_vlsr, $bad_fit_mask_22_vlsr, $bad_fit_mask_11_22_comb_vlsr);    #masks made from bad pixel lists 
    my ($at_least_one_good_11_fit, $at_least_one_good_22_fit,             # = 1 or 0 depending if there is at least one good fit in 
	$at_least_one_good_1122_fit);                                     #that transition
    my ($good_22_fit_at_11_max);                                          # is there a good 22 fit at max 11 pos?
    my ($good_11_fit_pixs, $good_22_fit_pixs, $good_1122_fit_pixs );      #list of pixel position with good fits for that transition
    my ($nelem_good_11_fit_pixs, $nelem_good_1122_fit_pixs);              #num pix with good 11 or 11+22 fits -> decide to map or not 
    my ($map_11_or_not, $map_1122_or_not);                                # if pix area > beam size then map that transition

    #variables which have had the bad fit mask applied
    # -> *_11bfm   = 11   bad fit mask applied
    # -> *_1122bfm = 1122 bad fit mask applied
    my ($hdr_11_11bfm);
    my ($nh3_11_amp_K_11bfm,   $nh3_11_amp_main_K_11bfm,   $nh3_11_fwhm_kms_11bfm,   $nh3_11_vel_kms_11bfm,   $nh3_11_resid_K_11bfm);
    my ($nh3_11_amp_K_1122bfm, $nh3_11_amp_main_K_1122bfm, $nh3_11_fwhm_kms_1122bfm, $nh3_11_vel_kms_1122bfm, $nh3_11_resid_K_1122bfm);
    my ($nh3_22_amp_K_1122bfm, $nh3_22_amp_main_K_1122bfm, $nh3_22_fwhm_kms_1122bfm, $nh3_22_vel_kms_1122bfm, $nh3_22_resid_K_1122bfm);

    #observed spectra, residuals and fit values for the pixel with the maximum reliable 11 fit amplitude
    my ($max_pix_nh3_11_amp_main_K_11bfm);  #getting pixel with the maximum reliable 11 fit amplitude
    my ($nh3_11_spect_K_at_max_pos_11bfm,   $nh3_11_resid_K_at_max_pos_11bfm);    #spectra and residual 
    my ($nh3_22_spect_K_at_max_pos_1122bfm, $nh3_22_resid_K_at_max_pos_1122bfm);  #spectra and residual 
    my ($nh3_11_amp_main_K_at_max_amp_11bfm,   $nh3_11_vel_kms_at_max_amp_11bfm,   $nh3_11_fwhm_kms_at_max_amp_11bfm);
    my ($nh3_22_amp_main_K_at_max_amp_1122bfm, $nh3_22_vel_kms_at_max_amp_1122bfm, $nh3_22_fwhm_kms_at_max_amp_1122bfm);

    #statistics for observed spectra (*samp*) and fit residuals (*ramp*) at the position of the
    #peak NH3(1,1) pixel with a reliable fit
    # - these stats are used to calculate the rms output to the latex table
    my ($mean_samp11,$prms_samp11,$median_samp11,$min_samp11,$max_samp11,$adev_samp11,$rms_samp11); #11 spect
    my ($mean_ramp11,$prms_ramp11,$median_ramp11,$min_ramp11,$max_ramp11,$adev_ramp11,$rms_ramp11); #11 resid
    my ($mean_samp22,$prms_samp22,$median_samp22,$min_samp22,$max_samp22,$adev_samp22,$rms_samp22); #22 spect
    my ($mean_ramp22,$prms_ramp22,$median_ramp22,$min_ramp22,$max_ramp22,$adev_ramp22,$rms_ramp22); #22 resid
    my ($sum_resid_11, $sum_spect_11, $sum_resid_22, $sum_spect_22); #used to see how much flux recovered

    #------------------------------
    #variables used in analysis
    my ($bad_analysis_pix, $analysis_mask);  #additional mask to deal with non-physica derived vals e.g. -ve tau

    #optical depth calculations
    my ($nh3_11_sat_rat_12_11bfm, $nh3_11_sat_rat_32_11bfm, #ratio of different combinations 
	$nh3_11_sat_rat_12_32_av_11bfm);                    #of nh3_11 satellite lines
    my ($nh3_11_sat_rat_12_32_av_11bfma);                   #spurious ratios removed (ie hyperine intensity anomalies)
    my ($nh3_11_tau_main_11bfma, $nh3_11_tau_tot_11bfma);   #nh3_11 main and total tau (after removing bad fits)
    my ($nh3_11_tau_main_1122bfma, $nh3_11_tau_tot_1122bfma); 
    my ($a_coeff, $b_coeff, $c_coeff, $d_coeff);            #coefficients used to calc tau from above ratios
    my ($good_tau_11_at_11_max);                            #is there a good fit for NH3_11 fit at the peak pixel?

    #column density calcs
    my ($A_JK_11, $beam_filling_factor, $T_bg);
    my ($col_dens_nh3_11_cm2,            $col_dens_nh3_tot_cm2); 
    my ($col_dens_nh3_11_cm2_at_max_amp, $col_dens_nh3_tot_cm2_at_max_amp); 

    #values derived at peak nh3_11 pixel that go into output table
    my ($good_analysis_at_max_amp,                          #is there a reliable value in analysis mask at peak 11? 
	$good_trot_at_max_amp,
	$good_tkin_at_max_amp,
	$good_col_dens_tot_at_max_amp,
	$good_fwhm_nt_at_max_amp,
	$nh3_11_tau_main_at_max_amp_11bfm,                  #nh3_11 main-line optical depth
	$t_rot_K_at_max_amp_1122bfm,                        #rotational temperature
	$t_kin_K_at_max_amp_1122bfm,                        #kinetic temperature
	$col_dens_tot_nh3_cm2_at_max_amp_1122bfm,           #column density
	$fwhm_nt_kms_at_max_amp_1122bfm);                   #non-thermal linewidth

    #pdl variables holding array of derived values
    my ($t_rot_K, $t_kin_K, $col_dens_tot_nh3_cm2, $fwhm_nontherm_kms, $fwhm_therm_kms);

    #------------------------------
    #variables for creating output spectra
    my ($spectra_ps_out_dev_11, $spectra_ps_out_name_11, $spectra_ps_out_dev_22, $spectra_ps_out_name_22);
    my ($spect_min_val_11, $res_min_val_11, $res_max_val_11, $spect_min_val_22, $res_min_val_22, $res_max_val_22);
    my ($spect_to_plot_11_vel_axis_kms, $spect_to_plot_22_vel_axis_kms);
    my (@spect_to_plot_11_vel_axis_kms_array );
    my (@spect_to_plot_22_vel_axis_kms_array);
    my (@nh3_11_vel_kms_at_max_amp_11bfm_array, @max_min_11, @vel_at_max_amp_22_array, @max_min_22, @sep_line_x, @sep_line_y);
 
    my ($spect_to_plot_11,       $spect_to_plot_22,       $resid_to_plot_11,             $resid_to_plot_22);
    my (@spect_to_plot_11_array, @spect_to_plot_22_array, @resid_to_plot_11_dcoff_array, @resid_to_plot_22_dcoff_array);
    my ($plot_resid_11, $plot_resid_22, $resid_to_plot_11_dcoff, $resid_to_plot_22_dcoff);
    my ($max_pos_data_11, $max_pos_data_22, $nh3_11_spect_at_max_pos_data_K, $nh3_22_spect_at_max_pos_data_K);
   
    #------------------------------
    #variable names for mapping

    #filenames for fits maps
    my ($amp_11_name_out,       $amp_22_name_out,         $vel_name_out, 
	$fwhm_name_out,         $resid_name_out);
    my ($amp_11_name_out_11bfm, $amp_22_name_out_1122bfm, $vel_name_out_11bfm, 
	$fwhm_name_out_11bfm,   $resid_name_out_11bfm);
    my ($tau_11_tot_name_out, $t_rot_name_out, $col_dens_tot_out);

    #miriad imaging variables
    my ($col_lookup, $amp_11_miriad_name, $amp_22_miriad_name, $vel_miriad_name);
    my ($fwhm_miriad_name, $tau_11_tot_miriad_name, $t_rot_miriad_name, $col_dens_tot_miriad_name);
    my ($output_name,$source_no_uscore, $label, $caption );
    my ($mean_amp,$prms_amp,$median_amp,$min_amp,$max_amp,$adev_amp,$rms_amp);
    my ($mean_amp22,$prms_amp22,$median_amp22,$min_amp22,$max_amp22,$adev_amp22,$rms_amp22);
    my ($mean_vel,$prms_vel,$median_vel,$min_vel,$max_vel,$adev_vel,$rms_vel);
    my ($mean_fwhm,$prms_fwhm,$median_fwhm,$min_fwhm,$max_fwhm,$adev_fwhm,$rms_fwhm);
    my ($mean_tau,$prms_tau,$median_tau,$min_tau,$max_tau,$adev_tau,$rms_tau);
    my ($mean_trot,$prms_trot,$median_trot,$min_trot,$max_trot,$adev_trot,$rms_trot);
    my ($mean_col,$prms_col,$median_col,$min_col,$max_col,$adev_col,$rms_col);
    #============================================================
    #READING INPUT FILE FROM CMD LINE AND SORTING FITTING PARS
    #
    #------------------------------
    #Get command line args
    $nfiles = scalar(@ARGV);
    if ($nfiles !=8) {
	die "\n\tIncorrect Input. Usage nh3_fit_sherpa_out.pl <SOURCE> 
                <INFILE1_BASE> <INFILE2_BASE> <HEADER_FILE_11> 
                <HEADER_FILE_22> <MASK_FILE_11> <MASK_FILE_22> region \n\n";
    }
    
    $source               = $ARGV[0];
    $nh3_11_amp_K_file    = $ARGV[1] . ".ampl.fits";
    $nh3_11_fwhm_kms_file = $ARGV[1] . ".fwhm.fits";
    $nh3_11_vel_kms_file  = $ARGV[1] . ".pos.fits";
    $nh3_11_resid_K_file  = $ARGV[1] . ".resid.fits";
    $nh3_22_amp_K_file    = $ARGV[2] . ".ampl.fits";
    $nh3_22_fwhm_kms_file = $ARGV[2] . ".fwhm.fits";
    $nh3_22_vel_kms_file  = $ARGV[2] . ".pos.fits";
    $nh3_22_resid_K_file  = $ARGV[2] . ".resid.fits";
    $header_file_11       = $ARGV[3];
    $header_file_22       = $ARGV[4];
    $duchamp_mask_file_11 = $ARGV[5];
    $duchamp_mask_file_22 = $ARGV[6];
    $region               = $ARGV[7];

    #------------------------------
    #Check files exist
    if ( (!(-e $nh3_11_amp_K_file ))   || (!(-e $nh3_22_amp_K_file ))  ||
	 (!(-e $nh3_11_fwhm_kms_file ))  || (!(-e $nh3_22_fwhm_kms_file ))  ||
	 (!(-e $nh3_11_vel_kms_file ))   || (!(-e $nh3_22_vel_kms_file ))   ||
	 (!(-e $nh3_11_resid_K_file )) || (!(-e $nh3_22_resid_K_file )) ||
	 (!(-e $header_file_11 ))     ||  (!(-e $header_file_22 ))    ||
	 (!(-e $duchamp_mask_file_11 ))     ||  (!(-e $duchamp_mask_file_22 ))
	) {
	die "Error reading files. Exiting. \n";
    }

    if ($region ne "none"){
	$region = "\"" . $region . "\"";
	print LOG "$region\n";
    }
    else{ 
	$region="";
	print LOG "No region being used.\n";
    }

    #------------------------------
    #reading in and checking files

    $hdr_11        = rfits($header_file_11,"hdrcopy");
    $hdr_22        = rfits($header_file_22,"hdrcopy");
    
    $nh3_11_amp_K    = rfits($nh3_11_amp_K_file,"hdrcopy");     
    $nh3_11_fwhm_kms = rfits($nh3_11_fwhm_kms_file,"hdrcopy");  
    $nh3_11_vel_kms  = rfits($nh3_11_vel_kms_file,"hdrcopy");                            
    $nh3_11_resid_K  = rfits($nh3_11_resid_K_file,"hdrcopy");        
    $nh3_11_duchamp_mask  = rfits($duchamp_mask_file_11,"hdrcopy"); 

    $nh3_22_amp_K    = rfits($nh3_22_amp_K_file,"hdrcopy");      
    $nh3_22_fwhm_kms = rfits($nh3_22_fwhm_kms_file,"hdrcopy"); 
    $nh3_22_vel_kms  = rfits($nh3_22_vel_kms_file,"hdrcopy");  
    $nh3_22_resid_K  = rfits($nh3_22_resid_K_file,"hdrcopy");    
    $nh3_22_duchamp_mask  = rfits($duchamp_mask_file_22,"hdrcopy"); 
    
    #converting to correct units
    $nh3_11_amp_K    *= 0.001; $nh3_22_amp_K    *= 0.001;       #converting from [mK] to [K]
    $nh3_11_fwhm_kms *= 0.001; $nh3_22_fwhm_kms *= 0.001;       #convert from [m/s] to [[km/s]
    $nh3_11_vel_kms  *= 0.001; $nh3_22_vel_kms  *= 0.001;       #convert from [m/s] to [[km/s]
    $nh3_11_vel_kms  += 19.6;                                   #change vlsr from outer sat to main comp

    $nh3_11_amp_main_K  = $nh3_11_amp_K->slice(":,:,(2)");  #using brackets selects only that plane of array
    $nh3_22_amp_main_K  = $nh3_22_amp_K->slice(":,:,(0)");  #now single gaussian not 5

    my ($tmp1,$tmp2);
    $tmp1=$nh3_11_amp_K->hdr->{NAXIS1}; $tmp2=$nh3_22_amp_K->hdr->{NAXIS1};
    print LOG "tmp1 = $tmp1  tmp2 = $tmp2 \n";
    $tmp1=$nh3_11_amp_K->hdr->{NAXIS2}; $tmp2=$nh3_22_amp_K->hdr->{NAXIS2};
    print LOG "tmp1 = $tmp1  tmp2 = $tmp2 \n";

    #checking file dimension are the same
    if( ($nh3_11_amp_K->hdr->{NAXIS1} != $nh3_22_amp_K->hdr->{NAXIS1}) || 
	($nh3_11_amp_K->hdr->{NAXIS2} != $nh3_22_amp_K->hdr->{NAXIS2})  ){
	print LOG "Error! NH3 files have different x,y dimensions.";
	exit;
    }
    
    #============================================================
    #APPLYING FIRST-PASS CUTS TO REMOVE BAD FITS AND CREATE RELIABLE FIT MASKS

    #------------------------------
    #getting rms -> using full residual cubes
    ($mean_ramp11_fullcube,$prms_ramp11_fullcube,$median_ramp11_fullcube,
	$min_ramp11_fullcube,$max_ramp11_fullcube,$adev_ramp11_fullcube,$rms_ramp11_fullcube) 
	= stats($nh3_11_resid_K);
    ($mean_ramp22_fullcube,$prms_ramp22_fullcube,$median_ramp22_fullcube,
	$min_ramp22_fullcube,$max_ramp22_fullcube,$adev_ramp22_fullcube,$rms_ramp22_fullcube) 
	= stats($nh3_22_resid_K);
    
    print LOG "fullcube $mean_ramp11_fullcube,$prms_ramp11_fullcube,$median_ramp11_fullcube,
               $min_ramp11_fullcube,$max_ramp11_fullcube,$adev_ramp11_fullcube,$rms_ramp11_fullcube \n";
    print LOG "fullcube $mean_ramp22_fullcube,$prms_ramp22_fullcube,$median_ramp22_fullcube,
               $min_ramp22_fullcube,$max_ramp22_fullcube,$adev_ramp22_fullcube,$rms_ramp22_fullcube \n";

    #outputting blind s/n maps
    $nh3_11_amp_main_K_sn_fullcube = $nh3_11_amp_main_K/$rms_ramp11_fullcube;
    $nh3_22_amp_main_K_sn_fullcube = $nh3_22_amp_main_K/$rms_ramp22_fullcube;

    $nh3_11_amp_main_K_sn_fullcube->wfits("nh3_11_amp_main_K_sn_fullcube.fits",-32);
    $nh3_22_amp_main_K_sn_fullcube->wfits("nh3_22_amp_main_K_sn_fullcube.fits",-32);

    #------------------------------
    #setting tolerances    
    $fwhm_tol_max_11_kms = 10.0;                         #no fits greater than 10kms
    $fwhm_tol_max_22_kms = 10.0;                         #no fits greater than 10kms
    $fwhm_tol_min_11_kms = $hdr_11->hdr->{CRPIX3}/1000;  #channel width 
    $fwhm_tol_min_22_kms = $hdr_22->hdr->{CRPIX3}/1000;  #channel width
    $vlsr_tol_11_kms     = 5.0;                          #Vlsr tol - after get peak spect
    $vlsr_tol_22_kms     = 5.0;                          #Vlsr tol - after get peak spect
    $amp_tol_max_11_K    = 5*max($hdr_11);               #dont accept fits >5*max in cube
    $amp_tol_max_22_K    = 5*max($hdr_22);               #dont accept fits >3*max in cube
    $amp_tol_min_11_K    = 5*$rms_ramp11_fullcube;       #min is 5*rms
    $amp_tol_min_22_K    = 5*$rms_ramp22_fullcube;       #min is 5*rms
    $min_pix_to_map      = 30;                           #need at least double the beam area 

    print LOG "fwhm_tol_max_11_kms=$fwhm_tol_max_11_kms, fwhm_tol_max_22_kms=$fwhm_tol_max_22_kms 
               fwhm_tol_min_11_kms=$fwhm_tol_min_11_kms fwhm_tol_min_22_kms=$fwhm_tol_min_22_kms \n";
    print LOG "amp_tol_max_11_K=$amp_tol_max_11_K amp_tol_max_22_K=$amp_tol_max_22_K 
               amp_tol_min_11_K=$amp_tol_min_11_K amp_tol_min_22_K=$amp_tol_min_22_K\n";

    #------------------------------
    #finding bad pixels based on the tolerances

    @nh3_11_fwhm_kms_dims   = dims($nh3_11_fwhm_kms);    #Previously these arrays had
    @nh3_11_amp_main_K_dims = dims($nh3_11_amp_main_K);  #different dimensions but 
    @nh3_22_fwhm_kms_dims   = dims($nh3_22_fwhm_kms);    #now they should be the same
    @nh3_22_amp_main_K_dims = dims($nh3_22_amp_main_K);

    print LOG "11 dims fwhm @nh3_11_fwhm_kms_dims  dims amp main @nh3_11_amp_main_K_dims \n";
    print LOG "22 dims fwhm @nh3_22_fwhm_kms_dims  dims amp main @nh3_22_amp_main_K_dims \n";
    
    print LOG "nh3_11_amp_main_K=$nh3_11_amp_main_K";
    print LOG "nh3_22_amp_main_K=$nh3_22_amp_main_K";

 

    $bad_fit_pix_11 = whichND( ($nh3_11_fwhm_kms < $fwhm_tol_min_11_kms) | 
				  ($nh3_11_fwhm_kms > $fwhm_tol_max_11_kms) | 
				  ($nh3_11_amp_main_K < $amp_tol_min_11_K) | 
				  ($nh3_11_amp_main_K > $amp_tol_max_11_K) #| 
			          #($bad_fit_mask_sherpa_11 < -9999)
                                   #($nh3_11_vel_kms < ($nh3_11_vel_at_max_amp_kms_11bfm - $vlsr_tol_11_kms) ) |
	                           #($nh3_11_vel_kms > ($nh3_11_vel_at_max_amp_kms_11bfm + $vlsr_tol_11_kms) )
                                  );

 #   $bad_fit_pix_22 = whichND( 
#			          ($bad_fit_mask_sherpa_22 < -100) 
#                                   
#	                         );
#    print LOG "bad_fit_pix_22 first = $bad_fit_pix_22 \n";

    $bad_fit_pix_22 = whichND( 
	                          ($nh3_22_fwhm_kms < $fwhm_tol_min_22_kms) | 
				  ($nh3_22_fwhm_kms > $fwhm_tol_max_22_kms) | 
				  ($nh3_22_amp_main_K < $amp_tol_min_22_K) | 
				  ($nh3_22_amp_main_K > $amp_tol_max_22_K) 
			          #($bad_fit_mask_sherpa_22 == 0) |
                                   #($nh3_22_vel_kms < ($nh3_11_vel_at_max_amp_kms_11bfm - $vlsr_tol) ) |
	                           #($nh3_22_vel_kms > ($nh3_11_vel_at_max_amp_kms_11bfm + $vlsr_tol) ) 
	                         );
    print LOG "bad_fit_pix_11 = $bad_fit_pix_11 \nbad_fit_pix_22 = $bad_fit_pix_22 \n";

    #------------------------------
    #creating bad pixel masks
    $bad_fit_mask_11 = $nh3_11_duchamp_mask;  #starting with emission region
    $bad_fit_mask_22 = $nh3_22_duchamp_mask;  #mask set by duchamp

    #outputting duchamp mask
    $nh3_11_duchamp_mask->wfits("duchamp_mask_11.fits",-32);
    $nh3_22_duchamp_mask->wfits("duchamp_mask_22.fits",-32);
    
    #outputting duchamp+sherpa mask
    $bad_fit_mask_sherpa_11 = $nh3_11_amp_K->slice(":,:,(2)")->sever;   #using amp as proxy to get bad sherpa fits 
    $bad_fit_mask_sherpa_22 = $nh3_22_amp_K->slice(":,:,(0)")->sever;   #using amp as proxy to get bad sherpa fits 

    badmask($bad_fit_mask_sherpa_11->inplace,0);
    badmask($bad_fit_mask_sherpa_22->inplace,0);
    
    $bad_fit_mask_sherpa_11->wfits("bad_fit_mask_sherpa_11f.fits",-32);
    $bad_fit_mask_sherpa_22->wfits("bad_fit_mask_sherpa_22f.fits",-32);

    print LOG "bad_fit_mask_sherpa_22=$bad_fit_mask_sherpa_22 \n";

    $good_fit_pix_sherpa_11 = whichND($bad_fit_mask_sherpa_11 != 0);
    $good_fit_pix_sherpa_22 = whichND($bad_fit_mask_sherpa_22 != 0);

    print LOG "good_fit_pix_sherpa_11=$good_fit_pix_sherpa_11 \ngood_fit_pix_sherpa_22=$good_fit_pix_sherpa_22 \n";

    $tmp = $bad_fit_mask_sherpa_11->indexND($good_fit_pix_sherpa_11);
    $tmp .=1;

    $tmp = $bad_fit_mask_sherpa_22->indexND($good_fit_pix_sherpa_22);
    $tmp .=1;

    print LOG "bad_fit_mask_sherpa_22=$bad_fit_mask_sherpa_22 \n";

    $bad_fit_mask_sherpa_11->wfits("bad_fit_mask_sherpa_11.fits",-32);
    $bad_fit_mask_sherpa_22->wfits("bad_fit_mask_sherpa_22.fits",-32);

    $bad_fit_mask_11 *= $bad_fit_mask_sherpa_11;
    $bad_fit_mask_22 *= $bad_fit_mask_sherpa_22;

    $bad_fit_mask_11->wfits("duchamp_sherpa_mask_11.fits",-32);
    $bad_fit_mask_22->wfits("duchamp_sherpa_mask_22.fits",-32);

    #print LOG "duchamp 11 mask \n $nh3_11_duchamp_mask \n duchamp 22 mask \n $nh3_22_duchamp_mask\n";

    $pix_change_11=0; $pix_change_22=0;
    unless($bad_fit_pix_11->isempty){	
	$tmp = $bad_fit_mask_11->indexND($bad_fit_pix_11);
	$tmp .=0;
	#$bad_fit_mask_11 = $bad_fit_mask_11->setbadif($bad_fit_mask_11==0);
	$bad_fit_mask_11->wfits("bad_fit_msk_11.fits",-32);
	#print LOG "changed 11 pix mask \n $bad_fit_mask_11 \n";
	$pix_change_11++;
    }
    unless($bad_fit_pix_22->isempty){
        $tmp = $bad_fit_mask_22->indexND($bad_fit_pix_22);
	$tmp .=0;
	#$bad_fit_mask_22 = $bad_fit_mask_22->setbadif($bad_fit_mask_22==0);
	$bad_fit_mask_22->wfits("bad_fit_msk_22.fits",-32);
	#print LOG "changed 22 pix mask \n $bad_fit_mask_22 \n";
	$pix_change_22++;
    }
    
    if($pix_change_11==0) {
	$bad_fit_mask_11->wfits("bad_fit_msk_11.fits",-32);
	print LOG "no pixels changed in 11 \n";
    }
    if($pix_change_22==0) {
	$bad_fit_mask_22->wfits("bad_fit_msk_22.fits",-32);
	print LOG "no pixels changed in 22 \n";
    }

    $bad_fit_mask_11_22_comb = $bad_fit_mask_11*$bad_fit_mask_22;
    $bad_fit_mask_11_22_comb->wfits("bad_fit_mask_11_22_comb.fits",-32);

    #============================================================
    #APPLYING BAD FIT MASKS TO VARIABLES

    $nh3_11_amp_K_11bfm        = $nh3_11_amp_K      * $bad_fit_mask_11;
    $nh3_11_amp_main_K_11bfm   = $nh3_11_amp_main_K * $bad_fit_mask_11;
    $nh3_11_fwhm_kms_11bfm     = $nh3_11_fwhm_kms   * $bad_fit_mask_11;
    $nh3_11_vel_kms_11bfm      = $nh3_11_vel_kms    * $bad_fit_mask_11;
    $nh3_11_resid_K_11bfm      = $nh3_11_resid_K    * $bad_fit_mask_11;

    $nh3_11_amp_K_1122bfm      = $nh3_11_amp_K      * $bad_fit_mask_11_22_comb;
    $nh3_11_amp_main_K_1122bfm = $nh3_11_amp_main_K * $bad_fit_mask_11_22_comb;
    $nh3_11_fwhm_kms_1122bfm   = $nh3_11_fwhm_kms   * $bad_fit_mask_11_22_comb;
    $nh3_11_vel_kms_1122bfm    = $nh3_11_vel_kms    * $bad_fit_mask_11_22_comb;
    $nh3_11_resid_K_1122bfm    = $nh3_11_resid_K    * $bad_fit_mask_11_22_comb;

    $nh3_22_amp_K_1122bfm      = $nh3_22_amp_K      * $bad_fit_mask_11_22_comb;
    $nh3_22_amp_main_K_1122bfm = $nh3_22_amp_main_K * $bad_fit_mask_11_22_comb;
    $nh3_22_fwhm_kms_1122bfm   = $nh3_22_fwhm_kms   * $bad_fit_mask_11_22_comb;
    $nh3_22_vel_kms_1122bfm    = $nh3_22_vel_kms    * $bad_fit_mask_11_22_comb;
    $nh3_22_resid_K_1122bfm    = $nh3_22_resid_K    * $bad_fit_mask_11_22_comb;
    
    #============================================================
    #DETERMINE WHICH ANALYSIS/IMAGING REGIME WE ARE IN FOR THIS SOURCE
    #BASED ON THE BAD FIT MASK PARAMETERS CALCULATED ABOVE

    #initialising as no good fits and no mapping
    $at_least_one_good_11_fit = $at_least_one_good_22_fit = $at_least_one_good_1122_fit = 0; 
    $good_22_fit_at_11_max = $map_11_or_not = $map_1122_or_not = 0;                             

    #checking for at least one good 11 fit?
    $good_11_fit_pixs = whichND($bad_fit_mask_11==1);                   #find any 11 good pixs
    unless($good_11_fit_pixs->isempty){$at_least_one_good_11_fit = 1;} 

    if($at_least_one_good_11_fit==1){

	#get pixel with maximum realiable 11 fit amplitude
	$max_pix_nh3_11_amp_main_K_11bfm = whichND( $nh3_11_amp_main_K_11bfm==max($nh3_11_amp_main_K_11bfm) );
	print LOG "max_pix_nh3_11_amp_main_K_11bfm=$max_pix_nh3_11_amp_main_K_11bfm \n";

	#enough pixels to map 11 emission? check if num good pixs 11 sig larger than beam area
	$nelem_good_11_fit_pixs = nelem($good_11_fit_pixs); #number of good fit elements (*2)
	                              #initialise with no mapping
	if( $nelem_good_11_fit_pixs/2 > $min_pix_to_map){   #enough good pixels to do mapping
	    $map_11_or_not = 1;
	    print LOG "map_11_or_not=$map_11_or_not; $nelem_good_11_fit_pixs/2 > $min_pix_to_map \n"
	}
	else{print LOG "map_11_or_not=$map_11_or_not; $nelem_good_11_fit_pixs/2 < $min_pix_to_map \n"}
 
	#now checking for good 22 
	$good_22_fit_pixs = whichND($bad_fit_mask_22==1);               #find any good 22 pixs
	unless($good_22_fit_pixs->isempty){$at_least_one_good_22_fit = 1;}

	if($at_least_one_good_22_fit == 1){ 

	    $good_1122_fit_pixs = whichND($bad_fit_mask_11_22_comb==1); #find any good pixs
	    unless($good_1122_fit_pixs->isempty){$at_least_one_good_1122_fit = 1;}
	    
	    #is there a good fit at max 11 pixel in 22 data? need () to get rid extra dims
	    $good_22_fit_at_11_max = $bad_fit_mask_22( ($max_pix_nh3_11_amp_main_K_11bfm->at(0,0)),
						       ($max_pix_nh3_11_amp_main_K_11bfm->at(1,0))  );
	    print LOG "good_22_fit_at_11_max=$good_22_fit_at_11_max \n";
	
	    #enough pixels to map 11+22 emission? check if num good pix 11+22 larger than beam area
	    $nelem_good_1122_fit_pixs = nelem($good_1122_fit_pixs); #number of good fit elements (*2)
	    if( $nelem_good_1122_fit_pixs/2 > $min_pix_to_map){     #enough good pixels to do mapping
		$map_1122_or_not = 1;
		print LOG "map_1122_or_not=$map_1122_or_not; $nelem_good_1122_fit_pixs/2 > $min_pix_to_map \n"
	    }
	    else{print LOG "map_1122_or_not=$map_1122_or_not; $nelem_good_1122_fit_pixs/2 < $min_pix_to_map \n"}
	}
    }
    
    print LOG "at_least_one_good_11_fit = $at_least_one_good_11_fit \n";
    print LOG "at_least_one_good_22_fit = $at_least_one_good_22_fit \n";
    print LOG "at_least_one_good_1122_fit = $at_least_one_good_1122_fit \n";

    #============================================================
    #NOW PROCEED DEPENDING ON WHICH PARTICULAR REGIME THIS SOURCE IS IN
    $good_analysis_at_max_amp = 0.0;
    #------------------------------------------------------------
    #REGIME 1: NO NH3(1,1) DETECTION
    if($at_least_one_good_11_fit==0){ print LOG "No NH3(1,1) detection.\n"; }

    elsif($at_least_one_good_11_fit==1){ #NH3(1,1) DETECTION 
	
	print LOG "At least one NH3(1,1) detection. \n";

	#------------------------------
	#finding location of peak NH3(1,1) emission and get spectra + fit resids 
	$max_pix_nh3_11_amp_main_K_11bfm  = whichND($nh3_11_amp_main_K_11bfm==max($nh3_11_amp_main_K_11bfm));
	$nh3_11_spect_K_at_max_pos_11bfm  = $hdr_11($max_pix_nh3_11_amp_main_K_11bfm->at(0,0), 
						    $max_pix_nh3_11_amp_main_K_11bfm->at(1,0),:)->sever;   
	$nh3_11_resid_K_at_max_pos_11bfm  = $nh3_11_resid_K($max_pix_nh3_11_amp_main_K_11bfm->at(0,0),
							    $max_pix_nh3_11_amp_main_K_11bfm->at(1,0),:)->sever;   

	$nh3_11_amp_main_K_at_max_amp_11bfm = $nh3_11_amp_main_K_11bfm->at($max_pix_nh3_11_amp_main_K_11bfm->at(0,0), 
								$max_pix_nh3_11_amp_main_K_11bfm->at(1,0));
	$nh3_11_vel_kms_at_max_amp_11bfm    = $nh3_11_vel_kms_11bfm->at($max_pix_nh3_11_amp_main_K_11bfm->at(0,0), 
								$max_pix_nh3_11_amp_main_K_11bfm->at(1,0));
	$nh3_11_fwhm_kms_at_max_amp_11bfm   = $nh3_11_fwhm_kms_11bfm->at($max_pix_nh3_11_amp_main_K_11bfm->at(0,0), 
								 $max_pix_nh3_11_amp_main_K_11bfm->at(1,0));

	print LOG "$max_pix_nh3_11_amp_main_K_11bfm \n";
	#print LOG "nh3_11_spect_K_at_max_pos_11bfm = $nh3_11_spect_K_at_max_pos_11bfm \n";
	#print LOG "nh3_11_resid_K_at_max_pos_11bfm = $nh3_11_resid_K_at_max_pos_11bfm \n";
	print LOG "nh3_11_vel_kms_at_max_amp_11bfm  = $nh3_11_vel_kms_at_max_amp_11bfm\n";
	print LOG "nh3_11_fwhm_kms_at_max_amp_11bfm = $nh3_11_fwhm_kms_at_max_amp_11bfm \n";
	
	#------------------------------
	#updating bad fit mask for fits with large vlsr offset from peak vlsr
	$bad_fit_mask_11_vlsr = $bad_fit_mask_11; 
 
	$bad_fit_pix_11_vlsr = whichND( ($nh3_11_vel_kms_11bfm   < ($nh3_11_vel_kms_at_max_amp_11bfm   - $vlsr_tol_11_kms) ) |
					($nh3_11_vel_kms_11bfm   > ($nh3_11_vel_kms_at_max_amp_11bfm   + $vlsr_tol_11_kms) ) );
	
	
	print LOG "$nh3_11_vel_kms_at_max_amp_11bfm   - $vlsr_tol_11_kms \n";
	print LOG "$nh3_11_vel_kms_at_max_amp_11bfm   + $vlsr_tol_11_kms\n";
	

	unless($bad_fit_pix_11_vlsr->isempty){	
	    $tmp = $bad_fit_mask_11_vlsr->indexND($bad_fit_pix_11_vlsr);
	    $tmp .=0;
	    $bad_fit_mask_11_vlsr->wfits("bad_fit_msk_11_vlsr.fits",-32);
	    $pix_change_11++;
	}
	
	#------------------------------
	#applying vlsr cut
		
	$nh3_11_amp_K_11bfm        = $nh3_11_amp_K_11bfm      * $bad_fit_mask_11_vlsr;
	$nh3_11_amp_main_K_11bfm   = $nh3_11_amp_main_K_11bfm * $bad_fit_mask_11_vlsr;
	$nh3_11_fwhm_kms_11bfm     = $nh3_11_fwhm_kms_11bfm   * $bad_fit_mask_11_vlsr;
	$nh3_11_vel_kms_11bfm      = $nh3_11_vel_kms_11bfm    * $bad_fit_mask_11_vlsr;
	$nh3_11_resid_K_11bfm      = $nh3_11_resid_K_11bfm    * $bad_fit_mask_11_vlsr;
	
	#------------------------------
        #getting statistics of nh3(1,1) spectra and resids
	($mean_samp11,$prms_samp11,$median_samp11,$min_samp11,$max_samp11,$adev_samp11,$rms_samp11) 
	    = stats($nh3_11_spect_K_at_max_pos_11bfm);
	($mean_ramp11,$prms_ramp11,$median_ramp11,$min_ramp11,$max_ramp11,$adev_ramp11,$rms_ramp11) 
	    = stats($nh3_11_resid_K_at_max_pos_11bfm);
	$sum_resid_11=sum($nh3_11_resid_K_at_max_pos_11bfm);
	$sum_spect_11=sum($nh3_11_spect_K_at_max_pos_11bfm);
	print LOG "peak s11 $mean_samp11,$prms_samp11,$median_samp11,$min_samp11,$max_samp11,
               $adev_samp11,$rms_samp11 \n";
	print LOG "peak r11 $mean_ramp11,$prms_ramp11,$median_ramp11,$min_ramp11,$max_ramp11,
               $adev_ramp11,$rms_ramp11 \n";
	open(STATS, ">>stats.txt") || die "Can't open file."; #outputting to stats file
	print STATS "$rms_ramp11 \t $sum_resid_11 \t $sum_spect_11 \t ";
	print STATS "$sum_resid_11 \t $sum_spect_11 \t\n";
	close(STATS);

	#------------------------------------------------------------
	#calculating tau 11 

	$analysis_mask = $bad_fit_mask_11;  #initialising analysis mask with 11 bad fit mask
	
	#getting intensity ratios of satellite lines
	$nh3_11_sat_rat_12_11bfm       = $nh3_11_amp_K_11bfm->slice(":,:,(1)")/$nh3_11_amp_K_11bfm->slice(":,:,(2)");
	$nh3_11_sat_rat_32_11bfm       = $nh3_11_amp_K_11bfm->slice(":,:,(3)")/$nh3_11_amp_K_11bfm->slice(":,:,(2)");
	$nh3_11_sat_rat_12_32_av_11bfm = ( $nh3_11_amp_K_11bfm->slice(":,:,(3)")+$nh3_11_amp_K_11bfm->slice(":,:,(1)") )/
	    (2*$nh3_11_amp_K_11bfm->slice(":,:,(2)"));
	
	#flagging analysis pixels if the ratio of inner/main is >1
	$bad_analysis_pix = whichND( $nh3_11_sat_rat_12_32_av_11bfm >=1);
	$tmp = $analysis_mask->indexND($bad_analysis_pix);
	$tmp .= 0;
	$nh3_11_sat_rat_12_32_av_11bfma = $nh3_11_sat_rat_12_32_av_11bfm * $analysis_mask;
	
	#using the following polynomial to get the tau from the brightness temp ratios of
	#inner and main satellite lines
	#f1(x)= a + b*x + c*(x**2) + d*(x**3)
	$a_coeff = -3.82407;         #+/- 0.01145      (0.2995%)
	$b_coeff = 18.4927;          #+/- 0.07087      (0.3832%)
	$c_coeff = -22.2117;         #+/- 0.1401       (0.6307%)
	$d_coeff = 17.3261;          #+/- 0.08896      (0.5135%)
	
	$nh3_11_tau_main_11bfma = $a_coeff + 
	    $b_coeff*$nh3_11_sat_rat_12_32_av_11bfma +
	    $c_coeff*$nh3_11_sat_rat_12_32_av_11bfma * $nh3_11_sat_rat_12_32_av_11bfma +
	    $d_coeff*$nh3_11_sat_rat_12_32_av_11bfma * $nh3_11_sat_rat_12_32_av_11bfma * $nh3_11_sat_rat_12_32_av_11bfma;
	
	#flagging analysis pixels if tau_main < 0
	$bad_analysis_pix = whichND( $nh3_11_tau_main_11bfma < 0);
	$tmp = $analysis_mask->indexND($bad_analysis_pix);
	$tmp .= 0;
	$nh3_11_tau_main_11bfma = $nh3_11_tau_main_11bfma * $analysis_mask;
	$nh3_11_tau_tot_11bfma  = $nh3_11_tau_main_11bfma / $frac_main_sat[0];
	
	#checking to see if nh3_11 peak pixel has a good tau value
	$good_tau_11_at_11_max = $analysis_mask( ($max_pix_nh3_11_amp_main_K_11bfm->at(0,0)),
	                                          ($max_pix_nh3_11_amp_main_K_11bfm->at(1,0))  );
	if($good_tau_11_at_11_max==1){	    
	    $nh3_11_tau_main_at_max_amp_11bfm = $nh3_11_tau_main_11bfma( ($max_pix_nh3_11_amp_main_K_11bfm->at(0,0)),
	                                          ($max_pix_nh3_11_amp_main_K_11bfm->at(1,0))  );
	    print LOG "good tau 11 val $nh3_11_tau_main_at_max_amp_11bfm\n";
	}
	else{ print LOG "bad tau 11 val \n";}
	print LOG "good_tau_11_at_11_max=$good_tau_11_at_11_max \n";

	#writing out fits files to check calcs 
	$nh3_11_sat_rat_12_32_av_11bfma->wfits("tau_11_ratio_main_sat.fits",-32); #ratio of main/sat
	$nh3_11_tau_main_11bfma->wfits("tau_11_main.fits",-32);
	$nh3_11_tau_tot_11bfma->wfits("tau_11_total.fits",-32);
	$analysis_mask->wfits("analysis_mask_tau11.fits",-32);

	#------------------------------------------------------------
	#REGIME 2A: AT LEAST ONE RELIABLE NH3(1,1) FIT BUT NO NH3(2,2) DETECTION
	if($good_22_fit_at_11_max==0){ print LOG "No good NH3(2,2) detection at peak 11 position.";}

	#------------------------------------------------------------
	#REGIME 2B: AT LEAST ONE RELIABLE NH3(1,1) FIT AND THERE IS A GOOD NH3(2,2) DETECTION
	#           AT THAT PIXEL TOO
	else{

	    print LOG "At least one good NH3(2,2) detection at the peak of the NH3(1,1) emission. \n";
	    $good_analysis_at_max_amp = 1.0;

	    #getting spectra and resids at position of peak NH3(1,1) emission
	    $nh3_22_spect_K_at_max_pos_1122bfm  = $hdr_22($max_pix_nh3_11_amp_main_K_11bfm->at(0,0), 
							$max_pix_nh3_11_amp_main_K_11bfm->at(1,0),:)->sever;   
	    $nh3_22_resid_K_at_max_pos_1122bfm  = $nh3_22_resid_K_1122bfm($max_pix_nh3_11_amp_main_K_11bfm->at(0,0),
								$max_pix_nh3_11_amp_main_K_11bfm->at(1,0),:)->sever;

	    $nh3_22_amp_main_K_at_max_amp_1122bfm = $nh3_22_amp_main_K_1122bfm->at($max_pix_nh3_11_amp_main_K_11bfm->at(0,0), 
								    $max_pix_nh3_11_amp_main_K_11bfm->at(1,0));
	    $nh3_22_vel_kms_at_max_amp_1122bfm  = $nh3_22_vel_kms_1122bfm->at($max_pix_nh3_11_amp_main_K_11bfm->at(0,0), 
								    $max_pix_nh3_11_amp_main_K_11bfm->at(1,0));
	    $nh3_22_fwhm_kms_at_max_amp_1122bfm = $nh3_22_fwhm_kms_1122bfm->at($max_pix_nh3_11_amp_main_K_11bfm->at(0,0), 
								     $max_pix_nh3_11_amp_main_K_11bfm->at(1,0));

	    #print LOG "nh3_22_spect_K_at_max_pos_1122bfm  = $nh3_22_spect_K_at_max_pos_1122bfm \n";
	    #print LOG "nh3_22_resid_K_at_max_pos_1122bfm  = $nh3_22_resid_K_at_max_pos_1122bfm \n";
	    print LOG "nh3_22_vel_kms_at_max_amp_1122bfm  = $nh3_22_vel_kms_at_max_amp_1122bfm \n";
	    print LOG "nh3_22_fwhm_kms_at_max_amp_1122bfm = $nh3_22_fwhm_kms_at_max_amp_1122bfm \n";


	    #------------------------------
	    #updating bad fit mask for fits with large vlsr offset from peak vlsr
	    $bad_fit_mask_22_vlsr = $bad_fit_mask_22; 
	    
	    $bad_fit_pix_22_vlsr = whichND( ($nh3_22_vel_kms_1122bfm   < ($nh3_22_vel_kms_at_max_amp_1122bfm   - $vlsr_tol_22_kms) ) |
					    ($nh3_22_vel_kms_1122bfm   > ($nh3_22_vel_kms_at_max_amp_1122bfm   + $vlsr_tol_22_kms) ) );
	    
	    
	    print LOG "$nh3_22_vel_kms_at_max_amp_1122bfm   - $vlsr_tol_22_kms \n";
	    print LOG "$nh3_22_vel_kms_at_max_amp_1122bfm   + $vlsr_tol_22_kms\n";
	    
	    
	    unless($bad_fit_pix_22_vlsr->isempty){	
		$tmp = $bad_fit_mask_22_vlsr->indexND($bad_fit_pix_22_vlsr);
		$tmp .=0;
		$bad_fit_mask_22_vlsr->wfits("bad_fit_msk_22_vlsr.fits",-32);
		$pix_change_22++;
	    }
	    
	    $bad_fit_mask_11_22_comb_vlsr = $bad_fit_mask_11_vlsr*$bad_fit_mask_22_vlsr;
	    $bad_fit_mask_11_22_comb_vlsr->wfits("bad_fit_mask_11_22_comb_vlsr.fits",-32);
	    
	    #------------------------------
	    #applying vlsr cut

	    $nh3_11_amp_K_1122bfm      = $nh3_11_amp_K_1122bfm      * $bad_fit_mask_11_22_comb_vlsr;
	    $nh3_11_amp_main_K_1122bfm = $nh3_11_amp_main_K_1122bfm * $bad_fit_mask_11_22_comb_vlsr;
	    $nh3_11_fwhm_kms_1122bfm   = $nh3_11_fwhm_kms_1122bfm   * $bad_fit_mask_11_22_comb_vlsr;
	    $nh3_11_vel_kms_1122bfm    = $nh3_11_vel_kms_1122bfm    * $bad_fit_mask_11_22_comb_vlsr;
	    $nh3_11_resid_K_1122bfm    = $nh3_11_resid_K_1122bfm    * $bad_fit_mask_11_22_comb_vlsr;
	    
	    $nh3_22_amp_K_1122bfm      = $nh3_22_amp_K_1122bfm      * $bad_fit_mask_11_22_comb_vlsr;
	    $nh3_22_amp_main_K_1122bfm = $nh3_22_amp_main_K_1122bfm * $bad_fit_mask_11_22_comb_vlsr;
	    $nh3_22_fwhm_kms_1122bfm   = $nh3_22_fwhm_kms_1122bfm   * $bad_fit_mask_11_22_comb_vlsr;
	    $nh3_22_vel_kms_1122bfm    = $nh3_22_vel_kms_1122bfm    * $bad_fit_mask_11_22_comb_vlsr;
	    $nh3_22_resid_K_1122bfm    = $nh3_22_resid_K_1122bfm    * $bad_fit_mask_11_22_comb_vlsr;
	    
	    #getting statistics of nh3(2,2) spectra and resids
	    ($mean_samp22,$prms_samp22,$median_samp22,$min_samp22,$max_samp22,$adev_samp22,$rms_samp22) 
		= stats($nh3_22_spect_K_at_max_pos_1122bfm);
	    ($mean_ramp22,$prms_ramp22,$median_ramp22,$min_ramp22,$max_ramp22,$adev_ramp22,$rms_ramp22) 
		= stats($nh3_22_resid_K_at_max_pos_1122bfm);
	    $sum_resid_22=sum($nh3_22_resid_K_at_max_pos_1122bfm);
	    $sum_spect_22=sum($nh3_22_spect_K_at_max_pos_1122bfm);
	    print LOG "peak s22 $mean_samp22,$prms_samp22,$median_samp22,$min_samp22,$max_samp22,
               $adev_samp22,$rms_samp22 \n";
	    print LOG "peak r22 $mean_ramp22,$prms_ramp22,$median_ramp22,$min_ramp22,$max_ramp22,
               $adev_ramp22,$rms_ramp22 \n";

	    #------------------------------------------------------------
	    #CALCULATIONS - following ungerechts et al 1986, A&A, 157, 207
	    
	    #initialising derived table vals to 1.0 for now
	    $t_rot_K_at_max_amp_1122bfm               = "1.0";
	    $t_kin_K_at_max_amp_1122bfm               = "1.0";
	    $fwhm_nt_kms_at_max_amp_1122bfm           = "1.0";
	    $col_dens_tot_nh3_cm2_at_max_amp_1122bfm  = "1.0";
	    
	    #updating analysis mask to include both reliable 11 and 22 fits
	    $analysis_mask = $analysis_mask * $bad_fit_mask_11_22_comb;

	    #------------------------------
	    #Trot - Ungerechts et al 1986, Eq A10, assuming: X11=X22, A11=A22 
	    $nh3_11_tau_tot_1122bfma = $nh3_11_tau_tot_11bfma * $analysis_mask; #updating tau11 mask
	    $t_rot_K = 1 - exp(-1*$nh3_11_tau_tot_1122bfma*$frac_main_sat[0]);
	    $t_rot_K = $t_rot_K * ($nh3_22_amp_main_K_1122bfm/$nh3_11_amp_main_K_1122bfm);
	    $t_rot_K = 1-$t_rot_K;
	    $t_rot_K = $t_rot_K->setbadif($t_rot_K<=0);             #flagging any pixels currently <= 0 
	    $analysis_mask = $analysis_mask->setbadif($t_rot_K<=0); #updating analysis mask
	    $t_rot_K = log($t_rot_K);;
	    $t_rot_K = -9 * $t_rot_K /(20*$nh3_11_tau_tot_1122bfma*$frac_main_sat[1]);
	    $t_rot_K = $t_rot_K->setbadif($t_rot_K<=0);             #flagging any pixels currently <= 0
	    $analysis_mask = $analysis_mask->setbadif($t_rot_K<=0); #updating analysis mask
	    $t_rot_K = log($t_rot_K);
	    $t_rot_K = -41/$t_rot_K;
	    $t_rot_K->wfits("trot.fits",-32);
	    #print LOG "nh3_11_tau_tot_1122bfma = $nh3_11_tau_tot_1122bfma \n";
	    #print LOG "t_rot_K = $t_rot_K \n";

	    #------------------------------
	    #nh3(1,1) column density  - Ungerechts et al 1986, Eq A13 
	    # - assuming term in square brackets = Tbg*bff
	    
	    $T_bg = 2.3;
	    $beam_filling_factor = 1.0;
	    $A_JK_11 = $nh3_11_amp_main_K_11bfm /(1 - exp(-1*$nh3_11_tau_main_11bfma));
	    
	    $col_dens_nh3_11_cm2 = 1.33773E13 *
		                         $nh3_11_tau_tot_11bfma *
					 $nh3_11_fwhm_kms_11bfm * # 1E4 *   #converting from [km/s] to [cm/s]
					 ($A_JK_11 + $beam_filling_factor*$T_bg);

	    $col_dens_nh3_11_cm2->wfits("col_dens_nh3_11_cm2.fits",-32);
	    $A_JK_11->wfits("AJK.fits", -32);

	    #------------------------------
	    #nh3 total column density  - Ungerechts et al 1986, Eq A15 
	    # - assuming term in square brackets = Tbg*bff

	    $col_dens_nh3_tot_cm2 = $col_dens_nh3_11_cm2 *
		                          ( 
					    (1/3)*exp(23.2/$t_rot_K) + 
					    1 +
					    (5/3)*exp(-41/$t_rot_K) +
					    (14/3)*exp(-100/$t_rot_K)
					  );
	
	    $col_dens_nh3_tot_cm2->wfits("col_dens_nh3_tot_cm2.fits", -32);
				  
	    #------------------------------
	    #tkin

	    $t_kin_K = $t_rot_K/(
		                  1 - 
		                  ($t_rot_K/42)*
		                  log(
			               1 + 1.1 * exp(-16/$t_rot_K)
                                     )
		               );
	    
	    $t_kin_K->wfits("t_kin_K.fits", -32);
  
	    #------------------------------
	    #delta v thermal and non-thermal

	    $fwhm_therm_kms = 0.23 * ($t_kin_K/20)**0.5 ;
	    $fwhm_therm_kms->wfits("fwhm_therm_kms.fits", -32);
	    
	    $fwhm_nontherm_kms = ($nh3_11_fwhm_kms_1122bfm*$nh3_11_fwhm_kms_1122bfm - $fwhm_therm_kms*$fwhm_therm_kms)**0.5;
	    $fwhm_nontherm_kms->wfits("fwhm_nontherm_kms.fits", -32);

	    #------------------------------
	    #determining whether peak pixel has a good value for the derived gas params
	    $good_trot_at_max_amp = $good_tkin_at_max_amp = 0.0;
	    $good_col_dens_tot_at_max_amp = $good_fwhm_nt_at_max_amp = 0;

	    #trot
	    $t_rot_K_at_max_amp_1122bfm  = $t_rot_K( ($max_pix_nh3_11_amp_main_K_11bfm->at(0,0)),
	                                          ($max_pix_nh3_11_amp_main_K_11bfm->at(1,0))  );
	    if( $t_rot_K_at_max_amp_1122bfm != 0 ){ $good_trot_at_max_amp = 1.0; }

	    #col_dens_11
	    $col_dens_nh3_11_cm2_at_max_amp = $col_dens_nh3_11_cm2( 
		                                            ($max_pix_nh3_11_amp_main_K_11bfm->at(0,0)),
	                                                    ($max_pix_nh3_11_amp_main_K_11bfm->at(1,0))  );
	    
	    #col_dens_tot
	    $col_dens_nh3_tot_cm2_at_max_amp = $col_dens_nh3_tot_cm2( 
		                                            ($max_pix_nh3_11_amp_main_K_11bfm->at(0,0)),
	                                                    ($max_pix_nh3_11_amp_main_K_11bfm->at(1,0))  );
	    if( $col_dens_nh3_tot_cm2_at_max_amp != 0){$good_col_dens_tot_at_max_amp = 1.0; }

	    #tkin
	    $t_kin_K_at_max_amp_1122bfm  = $t_kin_K( ($max_pix_nh3_11_amp_main_K_11bfm->at(0,0)),
	                                          ($max_pix_nh3_11_amp_main_K_11bfm->at(1,0))  );
	    if( $t_kin_K_at_max_amp_1122bfm != 0 ){ $good_tkin_at_max_amp = 1.0; }

	    #dvnt
	    $fwhm_nt_kms_at_max_amp_1122bfm = $fwhm_nontherm_kms( ($max_pix_nh3_11_amp_main_K_11bfm->at(0,0)),
	                                          ($max_pix_nh3_11_amp_main_K_11bfm->at(1,0))  );
	    if($fwhm_nt_kms_at_max_amp_1122bfm != 0){$good_fwhm_nt_at_max_amp = 1.0;}

	    print LOG "t_rot_K_at_max_amp_1122bfm            = $t_rot_K_at_max_amp_1122bfm \n";
	    print LOG "t_kin_K_at_max_amp_1122bfm            = $t_kin_K_at_max_amp_1122bfm \n";
	    print LOG "col_dens_nh3_11_cm2_at_max_amp  = $col_dens_nh3_11_cm2_at_max_amp \n";
	    print LOG "col_dens_nh3_tot_cm2_at_max_amp = $col_dens_nh3_tot_cm2_at_max_amp \n";
	}

    }
    
    #============================================================
    #NOW HAVE DERIVED VALUES, UPPER LIMITS OR NON DETECTIONS FOR EACH
    #PARAMETER, OUTPUT TO THE LATEXTABLE FILE

    #getting the source name
    $source_tex = $source; $source_tex =~ s/\_/ /g;
    $_ = $source_tex;
    if(/(G\S+)/){
	 $source_only = $1 ;
    }
    if(/G(\d+\.\d+)(\S+)/){
	 $l_deg = $1 ;
	 $b_deg = $2 ;
    }

    open (LATEXTABLE, ">>fits.tex") || die "Can't open latex file.";

    #16 columns 
    #print name and position no matter what
    print LATEXTABLE "$source_only  & $l_deg & $b_deg  ";

    #no 11 fits  %14 columns left
    if($at_least_one_good_11_fit==0){ #print non-detections for everything

	if(isbad($rms_ramp11_fullcube) && isbad($rms_ramp22_fullcube)){
	    #                  T11   V11   DV11 RMS11 Tau11 T22  V22   DV22   RMS22    Trot  N  Tkin DVnt 
	    print LATEXTABLE " & -   & -   & -  & -   & -   & -  & -   & -    & -      & -   & - & - & -  ";
	}
	elsif(isbad($rms_ramp11_fullcube)){
	    #                  T11   V11   DV11 RMS11 Tau11 T22  V22   DV22   RMS22    Trot  N  Tkin DVnt 
	    printf LATEXTABLE " & -   & -   & -  & -   & -   & -  & -   & -    & %4.2f    & -   & - & - & -  ",$rms_ramp22_fullcube;
	}
	elsif(isbad($rms_ramp22_fullcube)){
	    #                   T11   V11   DV11 RMS11     Tau11 T22  V22   DV22   RMS22  Trot  N  Tkin DVnt 
	    printf LATEXTABLE " & -   & -   & -  & %4.2f   & -   & -  & -   & -    & -    & -   & - & - & -  ",$rms_ramp11_fullcube;
	}
	else{
	    #                   T11   V11  DV11  RMS11   Tau11   T22   V22   DV22   RMS22    Trot  N  Tkin DVnt                 
	    printf LATEXTABLE " & -   & -   & -  & %4.2f   & -    & -    & -   & -   & %4.2f  & -   & - & - & -  ", $rms_ramp11_fullcube, $rms_ramp22_fullcube;
	}
    }

    else{ #good 11 fit
	#                 T11   V11     DV11     RMS11
	printf LATEXTABLE " & %4.2f & %4.2f & %4.2f & %4.2f   ", $nh3_11_amp_main_K_at_max_amp_11bfm,
                                                     $nh3_11_vel_kms_at_max_amp_11bfm,
                                                     $nh3_11_fwhm_kms_at_max_amp_11bfm,
	                                             $rms_ramp11_fullcube;
	
	if($good_tau_11_at_11_max==1){ #is there a good 11 tau value?
	    printf LATEXTABLE " & %4.2f   ", $nh3_11_tau_main_at_max_amp_11bfm;
	}
	else{print LATEXTABLE " & -   ";}  #bad 11 tau value at peak 11 pos
	
	


	if($good_22_fit_at_11_max==1){ #is there a good 22 fit at peak 11

	    #                     T22      V22    DV22     RMS22  
	    printf LATEXTABLE " & %4.2f & %4.2f & %4.2f &  %4.2f  ", $nh3_22_amp_main_K_at_max_amp_1122bfm,
                                                      $nh3_22_vel_kms_at_max_amp_1122bfm,
                                                      $nh3_22_fwhm_kms_at_max_amp_1122bfm,
                                                      $rms_ramp22_fullcube;

	    if($good_analysis_at_max_amp == 1){ #reliable value in analysis mask at peak 11? 

		$t_rot_K_at_max_amp_1122bfm->inplace->setnantobad;
		$col_dens_nh3_tot_cm2_at_max_amp->inplace->setnantobad;
		$t_kin_K_at_max_amp_1122bfm->inplace->setnantobad;
		$fwhm_nt_kms_at_max_amp_1122bfm->inplace->setnantobad;

		if($good_trot_at_max_amp !=0 && !(isbad($t_rot_K_at_max_amp_1122bfm)) ){
		    printf LATEXTABLE " & %4.0f ",$t_rot_K_at_max_amp_1122bfm;
		}
		else{                         print  LATEXTABLE " & - ";}

		if($good_col_dens_tot_at_max_amp !=0 && !(isbad($col_dens_nh3_tot_cm2_at_max_amp)) ){
		    
		    printf LATEXTABLE " & %4.2f ",($col_dens_nh3_tot_cm2_at_max_amp/1E15); 
		}
		else{                                 print  LATEXTABLE " & - ";}

		if($good_tkin_at_max_amp !=0 && !(isbad($t_kin_K_at_max_amp_1122bfm)) ){
		    printf LATEXTABLE " & %4.0f ",$t_kin_K_at_max_amp_1122bfm; 
		}
		else{                         print  LATEXTABLE " & - ";}

		if($good_fwhm_nt_at_max_amp !=0 && !(isbad($fwhm_nt_kms_at_max_amp_1122bfm)) ){
		    printf LATEXTABLE " & %4.2f ",$fwhm_nt_kms_at_max_amp_1122bfm; 
		}
		else{                            print  LATEXTABLE " & - ";}
		#                 Trot     N     Tkin     DVnt 
		#printf LATEXTABLE " & %4.2f & %4.2f & %4.2f & %4.2f ", $t_rot_K_at_max_amp_1122bfm,
                #                                      $t_kin_K_at_max_amp_1122bfm,
                #                                      $col_dens_tot_nh3_cm2_at_max_amp_1122bfm,
		#                                      $fwhm_nt_kms_at_max_amp_1122bfm;
	    }
	    else{
		#                    Trot   N  Tkin DVnt
		print LATEXTABLE " & -    & - & -   & - ";
	    }

	}
	else{ 
	    if(isbad($rms_ramp22_fullcube)){
		#                   T22   V22  DV22  RMS22 Trot N   Tkin DVnt
		printf LATEXTABLE " & -   & - & -   & -    & -  & - & -  & -";
	    }
	    else{
		#                   T22   V22  DV22  RMS22  Trot N   Tkin DVnt
		printf LATEXTABLE " & -   & - & -   & %4.2f & -  & - & -  & -", $rms_ramp22_fullcube;
	    }
	}
    }
    print LATEXTABLE " \\\\ \n";    
    close(LATEXTABLE);

    #============================================================
    #CREATING POSTSCRIPT FILES OF PEAK SPECTRA

    #creating ps device file name and output ps filename
    $spectra_ps_out_dev_11 = $source . "_nh3_11_spect.ps/ps";
    $spectra_ps_out_name_11 = $source . "_nh3_11_spect.ps";
    $spectra_ps_out_dev_22 = $source . "_nh3_22_spect.ps/ps";
    $spectra_ps_out_name_22 = $source . "_nh3_22_spect.ps";

    #choosing spectra and residual depending on reliability of detections
    if($at_least_one_good_11_fit == 0){

	$max_pos_data_11                = whichND($hdr_11==max($hdr_11));
	$max_pos_data_22                = whichND($hdr_22==max($hdr_22));
	$nh3_11_spect_at_max_pos_data_K = $hdr_11($max_pos_data_11->at(0,0), $max_pos_data_11->at(1,0),:)->sever; 
	$nh3_22_spect_at_max_pos_data_K = $hdr_22($max_pos_data_22->at(0,0), $max_pos_data_22->at(1,0),:)->sever; 

	$spect_to_plot_11 = $nh3_11_spect_at_max_pos_data_K;
	$resid_to_plot_11 = $nh3_11_spect_at_max_pos_data_K;
	$spect_to_plot_22 = $nh3_22_spect_at_max_pos_data_K;
	$resid_to_plot_22 = $nh3_22_spect_at_max_pos_data_K;
	$plot_resid_11    = 0;
	$plot_resid_22    = 0;

	#print LOG "spect_to_plot_11=$spect_to_plot_11 \nspect_to_plot_22=$spect_to_plot_22";
	#print LOG "resid_to_plot_11=$resid_to_plot_11 \nresid_to_plot_22=$resid_to_plot_22";
    }
    else{
	$plot_resid_11    = 1;
	$spect_to_plot_11 = $nh3_11_spect_K_at_max_pos_11bfm;
	$resid_to_plot_11 = $nh3_11_resid_K_at_max_pos_11bfm;
	
	if($good_22_fit_at_11_max==0){
	    $plot_resid_22    = 0;
	    $max_pos_data_22                = whichND($hdr_22==max($hdr_22));
	    $nh3_22_spect_at_max_pos_data_K = $hdr_22($max_pos_data_22->at(0,0), $max_pos_data_22->at(1,0),:)->sever;
	    
	    $spect_to_plot_22 = $nh3_22_spect_at_max_pos_data_K;
	    $resid_to_plot_22 = $nh3_22_spect_at_max_pos_data_K;
	}
	else{
	    $plot_resid_22    = 1;
	    $spect_to_plot_22 = $nh3_22_spect_K_at_max_pos_1122bfm;
	    $resid_to_plot_22 = $nh3_22_resid_K_at_max_pos_1122bfm;
	}
    }

    #getting axes ranges 
    $spect_min_val_11       = min($spect_to_plot_11); 
    $res_min_val_11         = min($resid_to_plot_11);
    $res_max_val_11         = max($resid_to_plot_11);
    $resid_to_plot_11_dcoff = $resid_to_plot_11 -
	$plot_resid_11 * ( sqrt($spect_min_val_11*$spect_min_val_11) +  sqrt($res_max_val_11*$res_max_val_11)   );

    $spect_min_val_22       = min($spect_to_plot_22); 
    $res_min_val_22         = min($resid_to_plot_22);
    $res_max_val_22         = max($resid_to_plot_22);
    $resid_to_plot_22_dcoff = $resid_to_plot_22 -
	$plot_resid_22 * ( sqrt($spect_min_val_22*$spect_min_val_22) +  sqrt($res_max_val_22*$res_max_val_22)   );

    #------------------------------
    #Calculating velocity of each pixel from header info
    $spect_to_plot_11_vel_axis_kms = 0.001 * (
	                           ( (sequence($hdr_11->hdr->{NAXIS3})
	                             - $hdr_11->hdr->{CRPIX3}) *
	                             $hdr_11->hdr->{CDELT3} ) +
	                           $hdr_11->hdr->{CRVAL3}
	                           );
    $spect_to_plot_22_vel_axis_kms = 0.001 * (
	                             ( (sequence($hdr_22->hdr->{NAXIS3})
	                               - $hdr_22->hdr->{CRPIX3}) *
	                               $hdr_22->hdr->{CDELT3} ) +
	                             $hdr_22->hdr->{CRVAL3}
	                             );

    #exit;

    #------------------------------
    #Creating ps files
    
    #converting pdl to array so can be read into pgplot
    @spect_to_plot_11_vel_axis_kms_array = list $spect_to_plot_11_vel_axis_kms;
    @spect_to_plot_11_array              = list $spect_to_plot_11;
    @resid_to_plot_11_dcoff_array        = list $resid_to_plot_11_dcoff;
    
    @spect_to_plot_22_vel_axis_kms_array = list $spect_to_plot_22_vel_axis_kms;
    @spect_to_plot_22_array              = list $spect_to_plot_22;
    @resid_to_plot_22_dcoff_array        = list $resid_to_plot_22_dcoff;	
    
    @max_min_11 = ( $spect_min_val_11-1.1*($res_max_val_11 - $res_min_val_11),max($hdr_11));
    @max_min_22 = ( $spect_min_val_22-1.1*($res_max_val_22 - $res_min_val_22),max($hdr_22));

    #
    if($at_least_one_good_11_fit == 0){
	print LOG "max_pos_data_11 = $max_pos_data_11  max_pos_data_22 = $max_pos_data_22 \n";
    }
    else{
	@nh3_11_vel_kms_at_max_amp_11bfm_array = ($nh3_11_vel_kms_at_max_amp_11bfm,$nh3_11_vel_kms_at_max_amp_11bfm);    
	@vel_at_max_amp_22_array = @nh3_11_vel_kms_at_max_amp_11bfm_array ;
	print LOG "max_pix_nh3_11_amp_main_K_11bfm = $max_pix_nh3_11_amp_main_K_11bfm \n";
	print LOG "nh3_11_vel_kms_at_max_amp_11bfm_array=@nh3_11_vel_kms_at_max_amp_11bfm_array \n";
	    
    };
    print LOG "plot_resid_11 = $plot_resid_11 plot_resid_22 = $plot_resid_22\n";
    print LOG "spect_min_val_11 = $spect_min_val_11: res_min_val_11 = $res_min_val_11: res_max_val_11=$res_max_val_11\n";
    print LOG "spect_min_val_22 = $spect_min_val_22: res_min_val_22 = $res_min_val_22: res_max_val_22=$res_max_val_22\n";
    #print LOG "spect_to_plot_11 = $spect_to_plot_11 \n";
    #print LOG "spect_to_plot_22 = $spect_to_plot_22 \n";
    #print LOG "resid_to_plot_11 = $resid_to_plot_11 \n";
    #print LOG "resid_to_plot_22 = $resid_to_plot_22 \n";
    #print LOG "resid_to_plot_11_dcoff = $resid_to_plot_11_dcoff \n";
    #print LOG "resid_to_plot_22_dcoff = $resid_to_plot_22_dcoff \n";
    #print LOG "spect_to_plot_11_vel_axis_kms = $spect_to_plot_11_vel_axis_kms \n";
    #print LOG "spect_to_plot_22_vel_axis_kms = $spect_to_plot_22_vel_axis_kms \n";
    #print LOG "spect_to_plot_11_array = @spect_to_plot_11_array \n";
    #print LOG "spect_to_plot_22_array = @spect_to_plot_22_array \n";
    #print LOG "resid_to_plot_11_dcoff_array = @resid_to_plot_11_dcoff_array \n";
    #print LOG "resid_to_plot_22_dcoff_array = @resid_to_plot_22_dcoff_array \n";
    #print LOG "max_min_11=@max_min_11 \n";
    #print LOG "max_min_22=@max_min_22 \n";

    #NH3 (1,1)
    pgopen($spectra_ps_out_dev_11);
    pgpap(10,1.0);
    pgsvp(0.13,0.99,0.01,0.99);
    pgwnad(0.0,1.0,0.0,0.7);
    pgslw(3.0);  pgscf(1); pgsch(1.3);
    pgswin($spect_to_plot_11_vel_axis_kms->at(0),$spect_to_plot_11_vel_axis_kms->at($hdr_11->hdr->{NAXIS3} - 1),
	   ( $spect_min_val_11-1.1*$plot_resid_11*($res_max_val_11 - $res_min_val_11)), max($hdr_11));
    pgtbox("bcnts",0.0,0.0,"vbncts",0.,0.);
    pgline($hdr_11->hdr->{NAXIS3},\@spect_to_plot_11_vel_axis_kms_array,\@spect_to_plot_11_array);
    pgline($hdr_11->hdr->{NAXIS3},\@spect_to_plot_11_vel_axis_kms_array,\@resid_to_plot_11_dcoff_array);
    pgsls(2);
    if($at_least_one_good_11_fit == 1){
	pgline(2, \@nh3_11_vel_kms_at_max_amp_11bfm_array, \@max_min_11);
    }
    pgsls(4);
    if($plot_resid_11){ #plotting line to separate spectra frum resid
	@sep_line_x = ($spect_to_plot_11_vel_axis_kms->at(0), $spect_to_plot_11_vel_axis_kms->at($hdr_11->hdr->{NAXIS3}-1));
	@sep_line_y = ($spect_min_val_11,$spect_min_val_11);
	pgline(2, \@sep_line_x, \@sep_line_y);
	print LOG "sep_line_x=@sep_line_x sep_line_y=@sep_line_y\n";
    }
    
    pgmtxt("L",3.0,0.5,0.5,"T\\dA\\u\\u*\\d [K]");
    pgmtxt("B",3.0,0.5,0.5,"V\\dLSR\\u [kms\\u-1\\d]");
    pgmtext("T",-1.4,0.05,0.0,"$source_only");
    pgmtext("T",-1.4,0.95,1.0,"NH\\d3\\u(1,1)");
    pgend();   

    #NH3 (2,2)
    pgopen($spectra_ps_out_dev_22);
    pgpap(10,1.0);
    pgsvp(0.13,0.99,0.01,0.99);
    pgwnad(0.0,1.0,0.0,0.7);
    pgslw(3.0);  pgscf(1); pgsch(1.3);
    pgswin($spect_to_plot_22_vel_axis_kms->at(0),$spect_to_plot_22_vel_axis_kms->at($hdr_22->hdr->{NAXIS3} - 1),
	   ( $spect_min_val_22-1.1*$plot_resid_22*($res_max_val_22 - $res_min_val_22)), max($hdr_22));
    pgtbox("bcnts",0.0,0.0,"vbncts",0.,0.);
    pgline($hdr_22->hdr->{NAXIS3},\@spect_to_plot_22_vel_axis_kms_array,\@spect_to_plot_22_array);
    pgline($hdr_22->hdr->{NAXIS3},\@spect_to_plot_22_vel_axis_kms_array,\@resid_to_plot_22_dcoff_array);
    pgsls(2);
    if($at_least_one_good_11_fit == 1){
	pgline(2, \@vel_at_max_amp_22_array, \@max_min_22);
    }
    pgsls(4);
    if($plot_resid_22){ #plotting line to separate spectra frum resid
	@sep_line_x = ($spect_to_plot_22_vel_axis_kms->at(0), $spect_to_plot_22_vel_axis_kms->at($hdr_22->hdr->{NAXIS3}-1));
	@sep_line_y = ($spect_min_val_22,$spect_min_val_22);
	pgline(2, \@sep_line_x, \@sep_line_y);
	print LOG "sep_line_x=@sep_line_x sep_line_y=@sep_line_y\n";
    }
    pgmtxt("L",3.0,0.5,0.5,"T\\dA\\u\\u*\\d [K]");
    pgmtxt("B",3.0,0.5,0.5,"V\\dLSR\\u [kms\\u-1\\d]");
    pgmtext("T",-1.4,0.05,0.0,"$source_only");
    pgmtext("T",-1.4,0.95,1.0,"NH\\d3\\u(2,2)");
    pgend();   
    

    #------------------------------
    #outputting appropriate info to latex file
    open (LATEXSPECTRA, ">>spectra_fig.tex") || die "Can't open file.";
    print LATEXSPECTRA "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$spectra_ps_out_name_11} &  \n";
    print LATEXSPECTRA "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$spectra_ps_out_name_22} \\\\  \n";
    close(LATEXSPECTRA);

    #============================================================
    #DECIDING WHETHER OR NOT TO MAP

    $col_lookup = -8;

    if($map_11_or_not == 1){ #yes, output nh3_11

	#------------------------------
	#output fits where bad fit mask applied
	$amp_11_name_out_11bfm   = $source . "_11_amp_main_11bfm.fits"; 
	$vel_name_out_11bfm      = $source . "_11_vel_main_11bfm.fits"; 
	$fwhm_name_out_11bfm     = $source . "_11_fwhm_main_11bfm.fits"; 
	$resid_name_out_11bfm    = $source . "_11_resid_wcs_11bfm.fits"; 
	
	badmask($nh3_11_amp_main_K_11bfm->inplace,0);
	badmask($nh3_11_vel_kms_11bfm->inplace,0);
	badmask($nh3_11_fwhm_kms_11bfm->inplace,0);
	badmask($nh3_11_resid_K_11bfm->inplace,0);

	$nh3_11_amp_main_K_11bfm   = $nh3_11_amp_main_K_11bfm->setbadif($nh3_11_amp_main_K_11bfm==0);
	$nh3_11_vel_kms_11bfm      = $nh3_11_vel_kms_11bfm->setbadif($nh3_11_vel_kms_11bfm==0);
	$nh3_11_fwhm_kms_11bfm     = $nh3_11_fwhm_kms_11bfm->setbadif($nh3_11_fwhm_kms_11bfm==0);
	$nh3_11_resid_K_11bfm      = $nh3_11_resid_K_11bfm->setbadif($nh3_11_resid_K_11bfm==0);

	#copying headers
	cp_hdr($nh3_11_amp_main_K_11bfm,   $hdr_11,"K");
	cp_hdr($nh3_11_vel_kms_11bfm,      $hdr_11,"km/s");    
	cp_hdr($nh3_11_fwhm_kms_11bfm,     $hdr_11,"km/s");
	cp_hdr($nh3_11_resid_K_11bfm,      $hdr_11,"K");
	
	$nh3_11_amp_main_K_11bfm->wfits($amp_11_name_out_11bfm,-32);		
	$nh3_11_vel_kms_11bfm->wfits($vel_name_out_11bfm,-32);
	$nh3_11_fwhm_kms_11bfm->wfits($fwhm_name_out_11bfm,-32);
	$nh3_11_resid_K_11bfm->wfits($resid_name_out_11bfm,-32);
	
	#------------------------------
	#getting limits
	($mean_amp,$prms_amp,$median_amp,$min_amp,$max_amp,$adev_amp,$rms_amp) 
	    = stats($nh3_11_amp_main_K_11bfm);
	($mean_vel,$prms_vel,$median_vel,$min_vel,$max_vel,$adev_vel,$rms_vel) 
	    = stats($nh3_11_vel_kms_11bfm);
	($mean_fwhm,$prms_fwhm,$median_fwhm,$min_fwhm,$max_fwhm,$adev_fwhm,$rms_fwhm) 
	    = stats($nh3_11_fwhm_kms_11bfm);
	
	print LOG "min_amp  = $min_amp   max_amp  = $max_amp \n";
	print LOG "min_vel  = $min_vel   max_vel  = $max_vel\n";
	print LOG "min_fwhm = $min_fwhm  max_fwhm = $max_fwhm\n";
	
	#------------------------------
	#miriad steps
	$amp_11_miriad_name=$amp_11_name_out_11bfm; 
	$amp_11_miriad_name =~s/fits/xy/;
	`fits in=$amp_11_name_out_11bfm out=$amp_11_miriad_name op=xyin`;
	print LOG "in=$amp_11_name_out_11bfm out=$amp_11_miriad_name op=xyin \n";
	
	$vel_miriad_name=$vel_name_out_11bfm; 
	$vel_miriad_name =~s/fits/xy/;
	`fits in=$vel_name_out_11bfm out=$vel_miriad_name op=xyin`;
	print LOG "in=$vel_name_out_11bfm out=$vel_miriad_name op=xyin \n";
	
	$fwhm_miriad_name=$fwhm_name_out_11bfm; 
	$fwhm_miriad_name =~s/fits/xy/;
	`fits in=$fwhm_name_out_11bfm out=$fwhm_miriad_name op=xyin`;
	print LOG "in=$fwhm_name_out_11bfm out=$fwhm_miriad_name op=xyin \n";
	
	#amp_11 on amp_11
	 `cgdisp in=$amp_11_miriad_name,$amp_11_miriad_name type=p,c slev=p,1 levs1=10,20,30,40,50,60,70,80,90 device=$source.amp_amp_11.ps/cps  beamtyp=b,l,2  options=wedge,blacklab labtyp=absdeg,absdeg csize=1.4,0,0,0 range=$min_amp,$max_amp,lin,$col_lookup cols1=1`;


	#vel_11 on amp_11 colour
	`cgdisp in=$vel_miriad_name,$amp_11_miriad_name type=p,c slev=p,1 levs1=10,20,30,40,50,60,70,80,90 device=$source.vel_amp.ps/cps  beamtyp=b,l,2  options=wedge,blacklab  labtyp=absdeg,absdeg csize=1.4,0,0,0 range=$min_vel,$max_vel,lin,$col_lookup cols1=1`;
	
	#fwhm_11 on amp_11
	`cgdisp in=$fwhm_miriad_name,$amp_11_miriad_name type=p,c slev=p,1 levs1=10,20,30,40,50,60,70,80,90 device=$source.fwhm_amp.ps/cps  beamtyp=b,l,2  options=wedge,blacklab  labtyp=absdeg,absdeg csize=1.4,0,0,0 range=$min_fwhm,$max_fwhm,lin,$col_lookup cols1=1`;
	
	#------------------------------
	#create latex file
	$output_name = $source . "output_fig.tex";
	$source_no_uscore = $source;
	$source_no_uscore =~ s/_/ /g;
	$source_no_uscore =~ s/img //g;
	$label="fig:" . $source . "_fig";
	$caption = "\\small{Figures for $source_only. NH\$_3\$(1,1) [top left] 
                   and NH\$_3\$(2,2) [top right] spectra with
                   residuals from the fit displayed underneath.
                   The left panel of the second row shows the integrated intensity 
                   emission of NH\$_3\$(1,1) 
                   overlayed with contours in 10\\% steps of the peak
                   value, from 90\\% down. The right panel of the second row
                   and left panel of the third tow show the FWHM
                   and V\$_{\\rm LSR}\$ of the NH\$_3\$(1,1) determined
                   by the fit. ";
                   
	open(OUTPUT, ">$output_name") || die "Can't open file.";
	print OUTPUT "\\documentclass[12pt,psfig]{mn2e} \n";
	print OUTPUT "\\usepackage{graphicx} \n";
	print OUTPUT "\\begin{document} \n";
	print OUTPUT "\\begin{figure*} \n";
	print OUTPUT "\\begin{center} \n";
	print OUTPUT "\\begin{tabular}{cc} \n";
	print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$spectra_ps_out_name_11} &  \n";
	print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$spectra_ps_out_name_22} \\\\  \n";
	print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source.amp_amp_11.ps} &  \n";
	print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source.fwhm_amp.ps} \\\\  \n";
	print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source.vel_amp.ps} &  \n";
	#print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source.amp_amp_22.ps} \\\\  \n";
#$source.amp_amp_22.ps

	if($map_1122_or_not==1){
	    
	    #applying bad masks
	    badmask($nh3_22_amp_main_K_1122bfm->inplace,0);
	    badmask($nh3_11_tau_tot_11bfma->inplace,0);
	    badmask($t_rot_K->inplace,0);
	    badmask($col_dens_nh3_tot_cm2->inplace,0);

	    $nh3_22_amp_main_K_1122bfm = $nh3_22_amp_main_K_1122bfm->setbadif($nh3_22_amp_main_K_1122bfm==0);	    
	    $nh3_11_tau_tot_11bfma=$nh3_11_tau_tot_11bfma->setbadif($nh3_11_tau_tot_11bfma==0);	    
	    $t_rot_K=$t_rot_K->setbadif($t_rot_K==0);	    
	    $col_dens_nh3_tot_cm2=$col_dens_nh3_tot_cm2->setbadif($col_dens_nh3_tot_cm2==0);

	    #------------------------------
	    #getting limits
	    ($mean_amp22,$prms_amp22,$median_amp22,$min_amp22,$max_amp22,$adev_amp22,$rms_amp22) 
		= stats($nh3_22_amp_main_K_1122bfm);
	    ($mean_tau,$prms_tau,$median_tau,$min_tau,$max_tau,$adev_tau,$rms_tau) 
		= stats($nh3_11_tau_tot_11bfma);
	    ($mean_trot,$prms_trot,$median_trot,$min_trot,$max_trot,$adev_trot,$rms_trot) 
		= stats($t_rot_K);
	    ($mean_col,$prms_col,$median_col,$min_col,$max_col,$adev_col,$rms_col) 
		= stats($col_dens_nh3_tot_cm2);
	
	    #print LOG "t_rot_K = $t_rot_K\n";
    
	    print LOG "min_amp22 = $min_amp22   max_amp22 = $max_amp22\n";
	    print LOG "min_tau   = $min_tau     max_tau   = $max_tau\n";
	    print LOG "min_trot  = $min_trot    max_trot  = $max_trot\n";
	    print LOG "min_col   = $min_col     max_col   = $max_col\n";

	    $tau_11_tot_name_out     = $source . "_11_tau_tot.fits"; 
	    $amp_22_name_out_1122bfm = $source . "_22_amp_main_1122bfm.fits"; 
	    $t_rot_name_out          = $source . "_t_rot.fits"; 
	    $col_dens_tot_out        = $source . "_col_dens_tot.fits";
	
	    cp_hdr($nh3_22_amp_main_K_1122bfm, $hdr_22,"K");
	    cp_hdr($nh3_11_tau_tot_11bfma,     $hdr_11,"Optical depth");
	    cp_hdr($t_rot_K,                   $hdr_11,"K");
	    cp_hdr($col_dens_nh3_tot_cm2,      $hdr_11,"cm-2");

	    $nh3_22_amp_main_K_1122bfm->wfits($amp_22_name_out_1122bfm,-32);
	    $nh3_11_tau_tot_11bfma->wfits($tau_11_tot_name_out,-32);
	    $t_rot_K->wfits($t_rot_name_out,-32);
	    $col_dens_nh3_tot_cm2->wfits($col_dens_tot_out,-32);

	    $amp_22_miriad_name=$amp_22_name_out_1122bfm; 
	    $amp_22_miriad_name =~s/fits/xy/;
	    `fits in=$amp_22_name_out_1122bfm out=$amp_22_miriad_name op=xyin`;
	    print LOG "in=$amp_22_name_out_1122bfm out=$amp_22_miriad_name op=xyin \n";
	    
	    $tau_11_tot_miriad_name=$tau_11_tot_name_out; 
	    $tau_11_tot_miriad_name =~s/fits/xy/;
	    `fits in=$tau_11_tot_name_out out=$tau_11_tot_miriad_name op=xyin`;
	    print LOG "in=$tau_11_tot_name_out out=$tau_11_tot_miriad_name \n";
	    
	    $t_rot_miriad_name=$t_rot_name_out; 
	    $t_rot_miriad_name =~s/fits/xy/;
	    `fits in=$t_rot_name_out out=$t_rot_miriad_name op=xyin`;
	    print LOG "in=$t_rot_name_out out=$t_rot_miriad_name \n";
	    
	    $col_dens_tot_miriad_name=$col_dens_tot_out; 
	    $col_dens_tot_miriad_name =~s/fits/xy/;
	    `fits in=$col_dens_tot_out out=$col_dens_tot_miriad_name op=xyin`;
	    print LOG "in=$col_dens_tot_out out=$col_dens_tot_miriad_name \n";

	    #amp_22 on amp_22
	    `cgdisp in=$amp_22_miriad_name,$amp_22_miriad_name type=p,c slev=p,1 levs1=10,20,30,40,50,60,70,80,90 device=$source.amp_amp_22.ps/cps  beamtyp=b,l,2  options=wedge,blacklab labtyp=absdeg,absdeg csize=1.4,0,0,0 range=$min_amp22,$max_amp22,lin,$col_lookup cols1=1 `;

	    #tau_11 on amp_11
	    `cgdisp in=$tau_11_tot_miriad_name,$amp_11_miriad_name type=p,c slev=p,1 levs1=10,20,30,40,50,60,70,80,90 device=$source.tau_amp.ps/cps  beamtyp=b,l,2  options=wedge,blacklab labtyp=absdeg,absdeg csize=1.4,0,0,0 range=$min_tau,$max_tau,lin,$col_lookup cols1=1`;

	    #t_rot on amp_11
	    `cgdisp in=$t_rot_miriad_name,$amp_11_miriad_name type=p,c slev=p,1 levs1=10,20,30,40,50,60,70,80,90 device=$source.t_rot_amp.ps/cps labtyp=absdeg,absdeg beamtyp=b,l,2  options=wedge,blacklab csize=1.4,0,0,0 range=$min_trot,$max_trot,lin,$col_lookup cols1=1`;
	    
	    #col dens on col dens
	     `cgdisp in=$col_dens_tot_miriad_name,$amp_11_miriad_name type=p,c slev=p,1 levs1=10,20,30,40,50,60,70,80,90 device=$source.col_dens_tot_amp.ps/cps labtyp=absdeg,absdeg beamtyp=b,l,2  options=wedge,blacklab csize=1.4,0,0,0 range=$min_col,$max_col,lin,$col_lookup cols1=1`;
	    
	    print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source.amp_amp_22.ps} \\\\  \n";
	    print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source.t_rot_amp.ps} &  \n";
	    print OUTPUT "\\includegraphics[height=6.5cm, angle=-90, trim=0 0 -5 0]{$source.col_dens_tot_amp.ps} \\\\  \n";
	    $caption = $caption . "The right panel on the third row shows the integrated intensity 
                   emission of NH\$_3\$(2,2) 
                   overlayed with contours in 10\\% steps of the peak
                   value, from 90\\% down. The bottom left and right panels show
                   the kinetic temperature and total NH\$_3\$ column density maps,
                   respectively.}"
	}
	else{$caption = $caption . "}";}

	print OUTPUT "\\end{tabular} \n";
	print OUTPUT "\\end{center} \n";
	print OUTPUT "\\caption{$caption}\n";
	print OUTPUT "\\label{$label} \n";
	print OUTPUT "\\end{figure*} \n";
	print OUTPUT "\\end{document} \n";
	close(OUTPUT);
	
	#------------------------------
	#running latex to create postscript

	`latex $output_name`; 
	$output_name=~s/tex/dvi/;   
	`dvips $output_name`;
	`rm *output_fig.aux *output_fig.dvi *output_fig.log *output_fig.tex`;

	#tidying up
	if($map_11_or_not == 1){`rm -r $amp_11_miriad_name $vel_miriad_name $fwhm_miriad_name`;}
	if($map_1122_or_not==1){`rm -r $tau_11_tot_miriad_name $t_rot_miriad_name $col_dens_tot_miriad_name $amp_22_miriad_name`;}

    }

    #============================================================
    #TIDYING UP

    if($map_11_or_not == 1){ 
	`rm $amp_11_name_out_11bfm $vel_name_out_11bfm $fwhm_name_out_11bfm $resid_name_out_11bfm `;
    }
    if($map_1122_or_not==1){`rm $tau_11_tot_name_out  $t_rot_name_out $col_dens_tot_out $amp_22_name_out_1122bfm`;}
    
    `rm $nh3_11_amp_K_file $nh3_11_fwhm_kms_file $nh3_11_vel_kms_file $nh3_11_resid_K_file`;
    `rm $nh3_22_amp_K_file $nh3_22_fwhm_kms_file $nh3_22_vel_kms_file $nh3_22_resid_K_file`;

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
