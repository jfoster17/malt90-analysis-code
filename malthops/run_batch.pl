#! /usr/bin/perl -w
#$Id $
#$Log$
#------------------------------------------------------------
#
#------------------------------------------------------------

#------------------------------
#Get command line args

#./run_batch.pl /mako4/b/malt90/malt90-pilot-data/hnc/ hnc 1
#./run_batch.pl /mako4/b/malt90/malt90-pilot-data/hcn/ hcn 1


$nfiles = scalar(@ARGV);
if ($nfiles !=3) {
    die "\n\tIncorrect Input. Usage run_batch.pl <DATA_DIR> <TRANSITION> <SMOOTHING>  \n\n"; 
}

$data_dir   = $ARGV[0];
$transition = $ARGV[1];
$pre_smooth = $ARGV[2];

#snl now only need to read in a single dir
#$data_dir = "/mako4/b/malt90/malt90-pilot-data/test/";

#use the following to clean up if have to abort code
#snl rm -r ./data/nh3_11/*_1.fits ./data/nh3_11/*.xy ./data/nh3_22/*_1.fits ./data/nh3_22/*.xy
#rm -r /mako4/b/malt90/malt90-pilot-data/test/*.xy /mako4/b/malt90/malt90-pilot-data/test/*_1.fits

#hanning smooth the data below if needed
#$pre_smooth = 1;

#snl keep all the latex stuff in for now as does no harm
#parameters for the latex table holding the fit values 
$columns_per_page = 40;        #number of columns in the table per page
$caption = "tmp";              #place holder text for table caption
$label = "tmp";                #place holder text for table label
$fig_rows_per_page = 4;        #number of rows of spectra per page in latex figure

#initialising column count to keep track of when to start new page in latex
$current_column_count = $current_spect_fig_cols_count = 1 ;

#============================================================
#READ IN AND PARSE ALL FILES

#read in all fits files from data dirs
@allfiles_dir=`ls $data_dir*fits`;

#------------------------------
#parsing files

$i=0;
foreach(@allfiles_dir){
    
    if(/$data_dir(.+)\.fits/){    #getting image files 
	$source[$i]=$1;
	print "$source[$i]\n";
	$i++;
    }
    else{print "didn't get one: $_ \n"; exit;} #bad filename
}

#snl need to keep track of total number of sources
$count1=$i;

#============================================================
#OPENING SOME TEXT FILES NEEDED TO HOLD SCRIPT OUTPUT

#------------------------------
##opening statistics file
if(-e "stats.txt"){
    `mv stats.txt stats_old.txt`;
}
open(STATS, ">stats.txt") || die "Can't open file.";
print STATS "#header\n";
print STATS "#rms_ramp_11 \t rms_ramp_22 \t sum_resid_11 \t sum_spect_11 \t ";
print STATS "sum_resid_22 \t sum_spect_22 \t \n";
close(STATS);

#------------------------------
#opening latex table file
if(-e "fits.tex"){
    `mv fits.tex fits_old.tex`;
}
open(LATEXTABLE, ">fits.tex") || die "Can't open file.";
print LATEXTABLE "\\documentclass[12pt,psfig]{mn2e} \n";
print LATEXTABLE "\\usepackage{graphicx} \n";
print LATEXTABLE "\\usepackage{rotating}\n";
print LATEXTABLE "\\begin{document} \n";
print LATEXTABLE "\\begin{sidewaystable*}[ht] \n";
print LATEXTABLE "\\begin{center} \n";
print LATEXTABLE "\\begin{tabular}{|c|c|c|c|c|} \n";  #5 cols
print LATEXTABLE "\\hline \\hline \n";
print LATEXTABLE "Source &  T\$\_{\\rm B}\$ & V\$\_{\\rm LSR}\$ & \$\\Delta\$V    & T\$\_{\\rm RMS}\$  \\\\ \n"; 
print LATEXTABLE "       & (K)              & (kms\$\^{-1}\$)   & (kms\$\^{-1}\$) & (K)                \\\\ \\hline \n";
close(LATEXTABLE);

#------------------------------------------------------------
#opening latex spectra file
if(-e "spectra_fig.tex"){
    `mv spectra_fig.tex spectra_fig_old.tex`;
}
open(LATEXSPECTRA, ">spectra_fig.tex") || die "Can't open file.";
print LATEXSPECTRA "\\documentclass[12pt,psfig]{mn2e} \n";
print LATEXSPECTRA "\\usepackage{graphicx} \n";
print LATEXSPECTRA "\\begin{document} \n";
print LATEXSPECTRA "\\begin{figure*} \n";
print LATEXSPECTRA "\\begin{center} \n";
print LATEXSPECTRA "\\begin{tabular}{cc} \n";
close(LATEXSPECTRA);

#------------------------------------------------------------
#GENERATING AND RUNNING BATCH FILE

#running through each source in the data dir
for($i=0;$i<$count1;$i++){
    
    #============================================================
    #RUNNING THE DATA FOR THIS SOURCE THROUGH THE FITTING PROCEDURE 

    #------------------------------
    #creating file that act as input file for fitting script
    open(OUTPUT, ">batch.txt") || die "Can't open file.";

    $src_base            = $source[$i] . "_sherpa";
    $src_name            = $data_dir . $source[$i] . ".fits";
    $src_mir_in_name = $src_name;  $src_mir_in_name=~s/fits/xy/;
    $src_mir_out_name = $data_dir . $source[$i] . "_" . $pre_smooth . ".xy";
    $src_mir_out_name_fits = $src_mir_out_name; $src_mir_out_name_fits=~s/xy/fits/;

    #snl the following lines set the sherpa fitting pars. I've hard coded them to fit
    #snl a single Gaussian with no mask file. We can change this later to fit multiple 
    #snl Gaussian's for transitions with hyperfine structure etc.
    #snl print OUTPUT "$src_base $src_name 5 True True ";
    #snl print OUTPUT "[-19600,-7450,0,7450,19600] $mask11_name None\n";
    #snl input pars: outname,fitsfile,ngauss,verbose,linkfwhm,linkpos,offsets,maskimg,plotname
    print OUTPUT    "$src_base $src_name 1 False False [0] None None \n";
    
    close(OUTPUT);

    #------------------------------
    #spectrally smoothing the data

    #print "fits in=$src_name out=$src_mir_in_name op=xyin \n";
    #print "hanning in=$src_mir_in_name out=$src_mir_out_name width=$pre_smooth \n";
    #print "fits in=$src_mir_out_name out=$src_mir_out_name_fits op=xyout \n";
    `fits in=$src_name out=$src_mir_in_name op=xyin`;
    `hanning in=$src_mir_in_name out=$src_mir_out_name width=$pre_smooth`;
    `fits in=$src_mir_out_name out=$src_mir_out_name_fits op=xyout`;

    #------------------------------
    #running "sherpa" fitting routine on the batch file created earlier
    # - this will output fits files with the fit results
    #   required by the analyis scripts in the next section

    `./fit_sherpa_batch_tl_hops.pl batch.txt`;   

    #exit;

    #============================================================
    #NOW HAVE FIT RESULTS, RUN THE ANALYSIS SCRIPTS

    #snl the following command line arguement was for the old version of 
    #snl the analysis fitting code to work with hops NH3 data.
    #snl new command line arguement belwo
    #this is the line that runs the analysis script "nh3_fit_sherpa_out_hops.pl"
    #snl ` ./nh3_fit_sherpa_out_hops.pl $nh3_11_source[$i] $src_11 $src_22 $src11_mir_out_name_fits $src22_mir_out_name_fits  $mask11_name $mask22_name none > tmp_nh3_fit.txt`;

    #snl new analysis file
    print "./fit_sherpa_out_malt.pl $source[$i] $src_base $src_name $transition\n";
    `./fit_sherpa_out_malt.pl $source[$i] $src_base $src_name $transition`;
    #snl ./fit_sherpa_out_malt.pl g305.36_NH3-c2h_MEAN g305.36_NH3-c2h_MEAN_sherpa /mako4/b/malt90/malt90-pilot-data/test/g305.36_NH3-c2h_MEAN.fits
    
    `rm -r $src_mir_in_name $src_mir_out_name $src_mir_out_name_fits`;

    print "$src_base sourcenum = $i: curr_col = $current_column_count: curr_spect_fig_cols = $current_spect_fig_cols_count\n";

    #------------------------------------------------------------
    #creating new page in latex table file if enough columns
    if($current_column_count >= $columns_per_page){
	
	#resetting column count
	$current_column_count = 0;

	#ending current table and starting newone 
	open(LATEXTABLE, ">>fits.tex") || die "Can't open file.";
	print LATEXTABLE "\\end{tabular} \n";
	print LATEXTABLE "\\end{center} \n";
        #print LATEXTABLE "\\caption{$caption}\n";  #will add caption manually later
	print LATEXTABLE "\\label{$label} \n";
        #print LATEXTABLE "\\end{table*} \n";       #commented this out to use sidewaystable instead
	print LATEXTABLE "\\end{sidewaystable*} \n\n";
	print LATEXTABLE "\\newpage \n \\clearpage \n\n";
	print LATEXTABLE "\\begin{sidewaystable*}[ht] \n";
	print LATEXTABLE "\\begin{center} \n";
	print LATEXTABLE "\\begin{tabular}{|c|c|c|c|c|} \n";  #5 cols
	print LATEXTABLE "\\hline \\hline \n";
	print LATEXTABLE "Source &  T\$\_{\\rm B}\$ & V\$\_{\\rm LSR}\$ & \$\\Delta\$V    & T\$\_{\\rm RMS}\$  \\\\ \n"; 
	print LATEXTABLE "       & (K)              & (kms\$\^{-1}\$)   & (kms\$\^{-1}\$) & (K)                \\\\ \\hline \n";
	close(LATEXTABLE);
    }
    
    #------------------------------------------------------------
    #creating new page in latex spectra file if enough columns
    if($current_spect_fig_cols_count >= $fig_rows_per_page){
	
	#resetting column count
	$current_spect_fig_cols_count = 0;

	#ending current table and starting newone 
	open(LATEXSPECTRA, ">>spectra_fig.tex") || die "Can't open file.";
	print LATEXSPECTRA "\\end{tabular} \n";
	print LATEXSPECTRA "\\end{center} \n";
        #print LATEXSPECTRA "\\caption{$caption}\n"; #will add caption manually later
	print LATEXSPECTRA "\\label{$label} \n";
        #print LATEXSPECTRA "\\end{table*} \n";
	print LATEXSPECTRA "\\end{figure*} \n\n";
	print LATEXSPECTRA "\\newpage \n \\clearpage \n\n";
	print LATEXSPECTRA "\\begin{figure*} \n";
	print LATEXSPECTRA "\\begin{center} \n";
	print LATEXSPECTRA "\\begin{tabular}{cc} \n";
	close(LATEXSPECTRA);
    }

    $current_column_count++;
    $current_spect_fig_cols_count++;
}

#------------------------------------------------------------
#closing latex table file and creating latex table ps
open(LATEXTABLE, ">>fits.tex") || die "Can't open file.";
print LATEXTABLE "\\hline \n";
print LATEXTABLE "\\end{tabular} \n";
print LATEXTABLE "\\end{center} \n";
#print LATEXTABLE "\\caption{$caption}\n";   #will add caption manually later
print LATEXTABLE "\\label{$label} \n";
#print LATEXTABLE "\\end{table*} \n";        #commented this out to use sidewaystable instead
print LATEXTABLE "\\end{sidewaystable*} \n";
print LATEXTABLE "\\end{document} \n";
close(LATEXTABLE);

#creating ps file
`latex fits.tex`; 
`dvips fits.dvi`;
`rm fits.aux fits.log fits.dvi `;

#------------------------------------------------------------
#closing latex spectra file and creating latex spectra ps
open(LATEXSPECTRA, ">>spectra_fig.tex") || die "Can't open file.";
print LATEXSPECTRA "\\end{tabular} \n";
print LATEXSPECTRA "\\end{center} \n";
#print LATEXSPECTRA "\\caption{$caption}\n";
print LATEXSPECTRA "\\label{$label} \n";
#print LATEXSPECTRA "\\end{table*} \n";
print LATEXSPECTRA "\\end{figure*} \n";
print LATEXSPECTRA "\\end{document} \n";
close(LATEXSPECTRA);

#creating ps file
`latex spectra_fig.tex`; 
`dvips spectra_fig.dvi`;
`rm spectra_fig.aux spectra_fig.log spectra_fig.dvi `;

#------------------------------------------------------------
#clearing up
if(!(-e "postscripts")){
	`mkdir postscripts`;
    }
    `mv *.ps postscripts`;

#getting rid of intermediate fits files used for debugging
#snl keep these for now
#`rm *fits *txt *tex`;
#`rm *fits *tex`;
`rm *fits`;

#use the following  line if want to make single "output.ps" from all mapped sources
#`ghostscript -dNOPAUSE -sDEVICE=pswrite -dBATCH -sOutputFile=output.ps postscripts/*output_fig.ps`;
#ghostscript -dNOPAUSE -sDEVICE=pswrite -dBATCH -sOutputFile=output.ps *output_fig.ps`

#some additional tidy up commands if required.
#rm -r *.xy *fits ; rm *ps ; rm postscripts/*ps
#rm *output_fig.aux *output_fig.dvi *output_fig.log *output_fig.tex

#rm -r /mako4/b/malt90/malt90-pilot-data/hcop/*xy /mako4/b/malt90/malt90-pilot-data/hcop/*_1.fits
#rm -r /mako4/b/malt90/malt90-pilot-data/test2/*xy /mako4/b/malt90/malt90-pilot-data/test2/*_1.fits

#============================================================
#Error messages
#./fit_sherpa_out_malt.pl g301_24um-hcop_MEAN g301_24um-hcop_MEAN_sherpa /mako4/b/malt90/malt90-pilot-data/hcop/g301_24um-hcop_MEAN.fits hcop
#Use of uninitialized value in numeric eq (==) at ./fit_sherpa_out_malt.pl line 465.
#Use of uninitialized value in concatenation (.) or string at ./fit_sherpa_out_malt.pl line 573.

