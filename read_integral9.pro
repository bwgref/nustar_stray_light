pro read_integral9, ra=ra, dec=dec, src_name, src_ra, src_dec, src_flag, slp_level=slp_level, nsrc=nsrc, Rmin=Rmin, Rmax=Rmax

  if(n_elements(Rmin) eq 0) then Rmin=1.0      
  if(n_elements(Rmax) eq 0) then Rmax=6.0      

  if(n_elements(slp_level) eq 0) then slp_level=4      

  if(slp_level eq 0) then begin
     fmin=0.0
  endif else if(slp_level eq 0) then begin
     fmin=3.0
  endif else if(slp_level eq 1) then begin
     fmin=5.0
  endif else if(slp_level eq 2) then begin
     fmin=10.0
  endif else if(slp_level eq 3) then begin
     fmin=50.0
  endif else if(slp_level eq 4) then begin
     fmin=100.0
  endif else begin
     print,'The max slp_level is 4'
     stop
  endelse
  if not (n_elements(flxmin) eq 0) then fmin=flxmin      

  cat=getenv('NUPLAN_AUXIL')+'/9year_integral_galsurvey.fits'
  
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
  if(selected eq 0) then begin 
     print,'No source selected '
     return
  endif else begin
     print,'Selected ',selected,' sources with 17-60 keV flux above ', fmin,' mCrab'     
     for t=0L, N_ELEMENTS(src_name)-1 do begin
        if not (src_flag[t] eq 1) then continue
        dist=sphdist(ra, dec, src_ra[t], src_dec[t], /DEGREES)
        
        if(dist lt Rmin or dist gt Rmax) then begin
           src_flag[t]=0
           ;print,'-',src_name[t],dist,' deg.'
           continue
        endif
        ;print,'+',src_name[t],src_flux[t],' mCrab',dist,' deg.',format='(2a,f8.4,a,f8.4,a)'
        print,src_name[t],' , ',src_flux[t],' , ',dist,format='(2a,f8.2,a,f8.4)'
     endfor
     print
  endelse
end
