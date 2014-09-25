pro read_combined_catalog, ra=ra, dec=dec, src_name, src_ra, src_dec, src_flux, src_flag, slp_level=slp_level, nsrc=nsrc, fmin=fmin

  ; fmin in mCrab
  if n_elements(fmin) eq 0 then fmin = 0. 

  Rmin=0.
  Rmax=5.0

  cat='auxil/9year_integral_BAT_70m_combined.fits'
  
  table=readfits(cat,header,EXTEN_NO=1,/SILENT)
  src_name=TBGET(header, table, 'NAME', /NOSCALE)
  src_ra=TBGET(header, table, 'RAJ2000', /NOSCALE)
  src_dec=TBGET(header, table, 'DEJ2000', /NOSCALE)
  src_flux=TBGET(header, table, 'FLUX1', /NOSCALE) ; 17-60 keV flux in mCrabs
  src_error=TBGET(header, table, 'E_FLUX1', /NOSCALE)

  src_flag=lonarr(N_ELEMENTS(src_name)) 
  src_flag(*)=0L
  index=where(src_flux gt fmin,selected)
  src_flag(index)=1L

  nsrc=selected
  openw, lun, /get_lun, 'straylight_sources.txt'
  if(selected eq 0) then begin 
     printf,lun,'No source selected '
     return
  endif else begin
     printf, lun, 'Selected ',selected,' sources with 17-60 keV flux above ', fmin,' mCrab'     
     for t=0L, N_ELEMENTS(src_name)-1 do begin
        if not (src_flag[t] eq 1) then continue
        dist=sphdist(ra, dec, src_ra[t], src_dec[t], /DEGREES)
        if(dist lt Rmin or dist gt Rmax) then begin
           src_flag[t]=0
           continue
        endif
        printf,lun,src_name[t],src_flux[t],' mCrab   ', dist, ' degrees away'
     endfor
     printf, lun, ""
  endelse
  goodones = where(src_flag eq 1, ngood)
  if ngood eq 0 then begin
     message, 'No stray light targets qualify!'
  endif else begin
     src_flag = src_flag[goodones]
     src_ra=src_ra[goodones]
     src_dec=src_dec[goodones]
     src_flux=src_flux[goodones]
     src_name=src_name[goodones]
  endelse

end
