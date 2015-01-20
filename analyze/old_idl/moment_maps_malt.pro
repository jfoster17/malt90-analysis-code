pro moment_maps_malt,directory
base_mom_maps = '/DATA/MALT_1/MALT90/data/moment_maps/'
in_dir = '/DATA/MALT_1/MALT90/data/gridzilla/'
cd,in_dir+directory,CURRENT=old_dir
files = file_search('*MEAN.fits')
for i=0,n_elements(files)-1 do begin
    fn = files[i]
    basename = file_basename(fn,'.fits')
    outdir = base_mom_maps+directory+'/'+basename+'_mom_maps'
    outname = basename
    if not file_test(outdir) then begin
        hd = headfits(fn)
                                ;CoordType = sxpar(hd, 'CTYPE1')
                                ;if CoordType eq 'RA---XXX' then begin
                                ;	print," **** Bad rotation in file **** "+fn
                                ;endif else begin
        print,"Processing... "+fn
        file_mkdir,outdir
        
        nplanes = sxpar(hd, 'NAXIS1')
        nrows   = sxpar(hd, 'NAXIS2')
        nchan   = sxpar(hd, 'NAXIS3')
        pad = 15
        vel_or = (findgen(nchan)+1-sxpar(hd, 'CRPIX3'))*sxpar(hd, 'CDELT3')+sxpar(hd, 'CRVAL3')
        
        mean = fltarr(nrows, nplanes)+!values.f_nan
        sd = fltarr(nrows, nplanes)+!values.f_nan
        errmn = fltarr(nrows, nplanes)+!values.f_nan
        errsd = fltarr(nrows, nplanes)+!values.f_nan
        skew = fltarr(nrows, nplanes)+!values.f_nan
        kurt = skew
        error = skew
        intint = skew
        npix = skew
        for z = 0, nplanes-1 do begin
            fits_in = readfits(fn,/silent)
            plane = reform(fits_in[z,*,*])
            ind = where(plane eq 0, ct)
            if ct gt 0 then plane[ind] = !values.f_nan
            goodspec = where(plane[*, 0] eq plane[*, 0], ct)
                                ;print,"Good Spectra:"
                                ;print,ct
            for y = 0, ct-1 do begin
                spec = plane[goodspec[y],*]
                if directory eq 'hc13ccn' or directory eq '13c34s' then begin
                   spec = spec[650:3500]
                   vel = vel_or[650:3500]
                   nchan = 2851
                endif else begin
					vel=vel_or
					spec[0:400] = 0
					spec[3695:4095] = 0
				endelse 
                djs_iterstat, spec, sigma = sigma
                spec2 = smooth(spec, 9, /edge_trun, /nan)
                spec2 = spec2[indgen(nchan/9)*9]
                mask = spec2 gt 1.5*sigma/3.
                mask = (shift(mask, 1)+shift(mask, -1))*mask gt 0
                mask[0] = 0b
                mask[n_elements(mask)-1] = 0b
                if total(mask) eq 0 then continue
                l = label_region(mask)
                h = histogram(l)
                null = max(h[1:*], thisone)
                mask = l eq thisone+1
                mask = congrid(mask, nchan)
                ind = where(mask)
                mask[((min(ind)-pad) > 0):((max(ind)+pad) < nchan-1)] = 1b
                ind = where(mask)
                mom = wt_moment(vel[ind], spec[ind], error = sigma+fltarr(n_elements(ind)))
                mean[goodspec[y], z]   = mom.mean
                sd[goodspec[y], z]     = mom.stdev
                errmn[goodspec[y], z]  = mom.errmn
                errsd[goodspec[y], z]  = mom.errsd
                skew[goodspec[y], z]   = total(spec[ind]*(vel[ind]-mom.mean)^3)/total(spec[ind])
                kurt[goodspec[y], z]   = total(spec[ind]*(vel[ind]-mom.mean)^4)/total(spec[ind])
                error[goodspec[y], z]  = sigma
                intint[goodspec[y], z] = total(spec[ind])
                npix[goodspec[y], z]   = n_elements(ind) 
            endfor
        endfor

        badind = where(errmn gt 1e3 or errsd gt 1e3, ct)
        if ct gt 0 then begin
            mean[badind] = !values.f_nan
            sd[badind] = !values.f_nan
        endif
        skew = skew/sd^3
        kurt = kurt/sd^4-3
        hdout = hd
        cd,outdir,CURRENT=tempdir
        
        error = reverse(rotate(error,1),1)
        mean = reverse(rotate(mean,1),1)
        sd = reverse(rotate(sd,1),1)
        errmn = reverse(rotate(errmn,1),1)
        errsd = reverse(rotate(errsd,1),1)
        kurt = reverse(rotate(kurt,1),1)
        skew = reverse(rotate(skew,1),1)
        intint = reverse(rotate(intint,1),1)
        
        sxaddpar, hdout, 'NAXIS', 2
        sxaddpar, hdout, 'NAXIS1', sxpar(hd, 'NAXIS1')
        sxaddpar, hdout, 'NAXIS2', sxpar(hd, 'NAXIS2')
        sxaddpar, hdout, 'CDELT1', sxpar(hd, 'CDELT1'), 'DEGREES'
        sxaddpar, hdout, 'CDELT2', sxpar(hd, 'CDELT2'), 'DEGREES'
        sxaddpar, hdout, 'CRPIX1', sxpar(hd, 'CRPIX1')
        sxaddpar, hdout, 'CRPIX2', sxpar(hd, 'CRPIX2')
        sxaddpar, hdout, 'CTYPE1', sxpar(hd, 'CTYPE1')
        sxaddpar, hdout, 'CTYPE2', sxpar(hd, 'CTYPE2')
        sxaddpar, hdout, 'CRVAL1', sxpar(hd, 'CRVAL1'), 'DEGREES'
        sxaddpar, hdout, 'CRVAL2', sxpar(hd, 'CRVAL2'), 'DEGREES'
        sxaddpar, hdout, 'BUNIT', 'KM/S'
        sxdelpar, hdout, 'CRVAL3'
        sxdelpar, hdout, 'CRPIX3'
        sxdelpar, hdout, 'CDELT3'
        sxdelpar, hdout, 'CTYPE3'
        writefits, outname+'.mom1.fits', float(mean)/1e3, hdout
        writefits, outname+'.mom2.fits', float(sd)/1e3, hdout
        writefits, outname+'.err1.fits', float(errmn)/1e3, hdout
        writefits, outname+'.err2.fits', float(errsd)/1e3, hdout
        sxaddpar, hdout, 'BUNIT', 'NONE'
        writefits, outname+'.mom3.fits', float(skew), hdout
        writefits, outname+'.mom4.fits', float(kurt), hdout

        sxaddpar, hdout, 'BUNIT', 'K.KM/S'
        intint = intint*sxpar(hd, 'CDELT3')/1e3
        writefits, outname+'.mom0.fits', float(intint), hdout
        writefits, outname+'.err0.fits', float(error)*sqrt(npix)*sxpar(hd, 'CDELT3')/1e3, hdout
        writefits, outname+'.emap.fits', float(error)*sxpar(hd, 'CDELT3')/1e3, hdout
        
        sxaddpar, hdout, 'BUNIT', 'SNR'
        writefits, outname+'.snr0.fits', float(intint)/float(error)*sqrt(npix)*sxpar(hd, 'CDELT3')/1e3,hdout
        
        cd,tempdir
                                ;endelse
    endif
endfor
cd,old_dir
return
end


