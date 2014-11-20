pro read_bat_5mcrab_catalog, ra=ra, dec=dec, src_name, src_ra, src_dec, src_flux, src_flag, slp_level=slp_level, nsrc=nsrc, fmin=fmin

  ; fmin in mCrab
  if n_elements(fmin) eq 0 then fmin = 0. 

  Rmin=0. ; Find everything from 0 to 5 degrees
  Rmax=5.0

  cat='auxil/BAT_70m_catalog_nustar_cut_5mcrab.fits'

  src = mrdfits(cat,1, /silent)


  src_name = src.name
  src_ra = src.raj2000
  src_dec = src.dej2000
  src_flux = src.flux1
  

  nsrc = n_elements(src) 


  ; All are already okay:
  src_flag=lonarr(nsrc)
;  index=where(src_flux gt fmin,selected)
;  src_flag(index)=1L
  src_flag(*) = 1L
;  nsrc=selected

  openw, lun, /get_lun, 'straylight_sources.txt'
  
;  printf, lun, 'Catalog has ' sources with 3-30 keV flux above 5 mCrab'     
  for t=0L, nsrc-1 do begin
     dist=sphdist(ra, dec, src_ra[t], src_dec[t], /DEGREES)
     if(dist lt Rmin or dist gt Rmax) then begin
        src_flag[t]=0
        continue
     endif
     outstring = src_name[t] + '  '+string(src_flux[t], format = '(d8.3)')+ ' mCrab   '+string(dist, format = '(d8.3)')+' away'
     printf, lun, outstring
  endfor
  printf, lun, ""
  
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
