
; e.g.
;root='/Data/NuSTAR/data/mini-survey'
;obsid_array=['40010001002', '40010002001', '40010003001', '40010004001', '40010005001', '40010006001']

pro make_sky_pdf,root=root,obsid=obsid,silent=silent, e1=e1,e2=e2, grade_max=grade_max, evtpath=evtpath, PA=PA, key=key, sky=sky, hdr=hdr, delta=delta, igne1=igne1, igne2=igne2, center_ra=center_ra, center_dec=center_dec

  skyw=1000L
  det1w=360
  hotpix_max=1000
  cha1=(e1-1.6)/0.04
  cha2=(e2-1.6)/0.04

  if(n_elements(igne1) eq 0) then begin
     ign1=(e2-1.6)/0.04
     ign2=(e1-1.6)/0.04
  endif else begin
     ign1=(igne1-1.6)/0.04
     ign2=(igne2-1.6)/0.04
  endelse

  ab=['A','B']
  npix=64
  hotpix_raw=lonarr(npix,npix,2)

  if(n_elements(evtpath) eq 0) then evtpath='event_cl'
  if(n_elements(key) eq 0) then key='nuplan'

  if(n_elements(ignore) eq 0) then ignore=[e2,e1]

  if(n_elements(root) eq 0) then begin
     print,"Give me root directory! root='path to data'" 
     return
  endif
  if(n_elements(obsid) eq 0) then begin
     print,"Give me obsid! e.g. obsid='40010001002'"
     return
  endif


  sky=DBLARR(skyw,skyw,2)

  files_evt = root+'/'+obsid+['/'+evtpath+'/nu'+obsid+'A01_cl_pdf.evt','/'+evtpath+'/nu'+obsid+'B01_cl_pdf.evt']
  ;get_raw_area,root=root,obsid=obsid,/silent,eventsdir=evtpath,area=area,mask_det1=mask_det1,mask_raw=mask_raw
     
  for t=0,1 do begin
     if(~file_test(files_evt[t])) then continue
     print,'Read ',files_evt[t]
        
     table=readfits(files_evt[t],header,EXTEN_NO=1,/SILENT)
    
     ;grade=TBGET( header, table, 'GRADE', /NOSCALE)
     chan=TBGET( header, table, 'PI', /NOSCALE)
     det_id=TBGET( header, table, 'DET_ID', /NOSCALE)
     det1x=TBGET( header, table, 'DET1X', /NOSCALE)-1
     det1y=TBGET( header, table, 'DET1Y', /NOSCALE)-1
     xsky=TBGET( header, table, 'X', /NOSCALE)
     ysky=TBGET( header, table, 'Y', /NOSCALE)
     rawx=TBGET( header, table, 'RAWX', /NOSCALE)
     rawy=TBGET( header, table, 'RAWY', /NOSCALE)
     weight=TBGET( header, table, 'WEIGHT', /NOSCALE)
     nx = SXPAR(Header,'TLMAX13',COUNT=A)
     ny = SXPAR(Header,'TLMAX14',COUNT=A)
     live = SXPAR(Header,'LIVETIME',COUNT=A)
     crval1 = SXPAR(Header,'TCRVL13',COUNT=A)
     crval2 = SXPAR(Header,'TCRVL14',COUNT=A)
     crpix1 = SXPAR(Header,'TCRPX13',COUNT=A)
     crpix2 = SXPAR(Header,'TCRPX14',COUNT=A)
     cdelt1 = SXPAR(Header,'TCDLT13',COUNT=A)
     cdelt2 = SXPAR(Header,'TCDLT14',COUNT=A)
     ctype1 = SXPAR(Header,'TCTYP13',COUNT=A)
     ctype2 = SXPAR(Header,'TCTYP14',COUNT=A)
     PA = SXPAR(Header,'PA_PNT',COUNT=A)

     delta=[cdelt1, cdelt2]
     make_astr,astr, DELTA = [cdelt1, cdelt2], $
               CRPIX = [crpix1,crpix2], CRVAL = [crval1, crval2], CTYPE=[ctype1,ctype2],$
               RADECSYS = 'FK5', EQUINOX = 2000.0
     center_ra=crval1
     center_dec=crval2
     mkhdr, hdr, 3, [skyw,skyw]
     putast, hdr, astr, CD_TYPE=2
     SXADDPAR, hdr, 'LIVETIME', live, ' Exposure, s'
     SXADDPAR, hdr, 'EMIN', e1, ' EMIN, keV'
     SXADDPAR, hdr, 'EMAX', e2, ' EMIN, keV'
     SXADDPAR, hdr, 'OBSID', obsid, ''
     SXADDPAR, hdr, 'MODULE', ab[t], ' Module'

     for i=0L,N_ELEMENTS(xsky)-1 do begin
        if (chan[i] lt cha1 or chan[i] gt cha2) then continue
        if (chan[i] gt ign1 and chan[i] lt ign2) then continue
        if(xsky[i] ge 0 and ysky[i] ge 0) then begin
           sky[xsky[i],ysky[i],t]+=weight[i]           
        endif
     endfor

     writefits,obsid+ab[t]+'.'+key+'.sky_pdf.fits',sky[*,*,t],hdr

  endfor
end
