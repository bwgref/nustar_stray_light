pro read_3to30_catalog, ra=ra, dec=dec, src_name, src_ra, src_dec, src_flux, src_flag, int_src, slp_level=slp_level, nsrc=nsrc, fmin=fmin

  ; fmin in mCrab
  if n_elements(fmin) eq 0 then fmin = 0. 

  Rmin=0. ; Find everything from 0 to 5 degrees
  Rmax=5.0

  cat='auxil/bat_integral_3to30keV.fits'

  src = mrdfits(cat,1, /silent)



  src_name = src.name
  src_ra = src.ra
  src_dec = src.dec
  src_flux = src.flux
  int_src = src.int_src ; is this an INTEGRAL only source? (1 = yes, 0 = no)
  
  ;; table=readfits(cat,header,EXTEN_NO=1,/SILENT)
  ;; src_name=TBGET(header, table, 'NAME', /NOSCALE)
  ;; src_ra=TBGET(header, table, 'RAJ2000', /NOSCALE)
  ;; src_dec=TBGET(header, table, 'DEJ2000', /NOSCALE)
  ;; src_flux=TBGET(header, table, 'FLUX1', /NOSCALE) ; 17-60 keV flux in mCrabs
  ;; src_error=TBGET(header, table, 'E_FLUX1', /NOSCALE)

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
        if int_src[t] eq 1.0 then warn = ' *** flux in INTEGRAL 17-60 keV band!' else warn = ''
        outstring = src_name[t] + '  '+string(src_flux[t], format = '(d8.3)')+ ' mCrab   '+string(dist, format = '(d8.3)')+' away'+warn
        printf, lun, outstring
;        printf,lun,src_name[t],src_flux[t],' mCrab   ', dist, ' degrees away'+warn
     endfor
     printf, lun, ""
  endelse
  close, lun
  free_lun, lun

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
