
; e.g.
;root='/Data/NuSTAR/data/mini-survey'
;obsid_array=['40010001002', '40010002001', '40010003001', '40010004001', '40010005001', '40010006001']

pro make_sky_oa,root=root,obsid=obsid,silent=silent, e1=e1,e2=e2, grade_max=grade_max, evtpath=evtpath, PA=PA, key=key, $
                optsky=optsky, sky=sky, hdr=hdr, delta=delta, igne1=igne1, igne2=igne2, oax_mm=oax_mm, oay_mm=oay_mm, $
                ra_obj=ra_obj,dec_obj=dec_obj

  skyw=1000L
  det1w=360
  hotpix_max=1000
  cha1=(e1-1.6)/0.04
  cha2=(e2-1.6)/0.04

  npix=64
  optraw=lonarr(npix,npix,2)  


  if(n_elements(igne1) eq 0) then begin
     ign1=(e2-1.6)/0.04
     ign2=(e1-1.6)/0.04
  endif else begin
     ign1=(igne1-1.6)/0.04
     ign2=(igne2-1.6)/0.04
  endelse

  oax_mm=DBLARR(2)
  oay_mm=DBLARR(2)


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


  sky=LONARR(skyw,skyw,2)
  optsky=LONARR(skyw,skyw,2)
  optdet=LONARR(det1w,det1w,2)

  files_evt = root+'/'+obsid+['/'+evtpath+'/nu'+obsid+'A01_cl.evt','/'+evtpath+'/nu'+obsid+'B01_cl.evt']
  files_opt = root+'/'+obsid+['/'+evtpath+'/nu'+obsid+'A_oa.fits','/'+evtpath+'/nu'+obsid+'B_oa.fits']
  files_mst = root+'/'+obsid+'/'+evtpath+'/nu'+obsid+'_mast.fits'
  files_att = root+'/'+obsid+'/'+evtpath+'/nu'+obsid+'_att.fits'
  ;get_raw_area,root=root,obsid=obsid,/silent,eventsdir=evtpath,area=area,mask_det1=mask_det1,mask_raw=mask_raw
     
  print,'***'
  for t=0,1 do begin
     if(~file_test(files_evt[t])) then begin 
        print,'? ',files_evt[t]
        continue 
     endif
     if(~file_test(files_opt[t])) then continue
     if(~file_test(files_att)) then continue
     if(~file_test(files_mst)) then continue
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
     RA_OBJ = SXPAR(Header,'RA_OBJ',COUNT=A)
     DEC_OBJ = SXPAR(Header,'DEC_OBJ',COUNT=A)
     print,'RA_OBJ,DEC_OBJ: ',RA_OBJ,DEC_OBJ

     delta=[cdelt1, cdelt2]
     make_astr,astr, DELTA = [cdelt1, cdelt2], CRPIX = [crpix1,crpix2], CRVAL = [crval1, crval2], CTYPE=[ctype1,ctype2]
               ;RADECSYS = 'FK5', EQUINOX = 2000.0
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
           sky[xsky[i],ysky[i],t]++
        endif
     endfor

     writefits,obsid+ab[t]+'.'+key+'.sky.fits',sky[*,*,t],hdr

     fpd = lonarr(npix/2.,npix/2.,4)
     for i=0L,N_ELEMENTS(rawx)-1 do begin
        if(rawx[i] ge 0 and rawy[i] ge 0) then begin
           fpd[rawx[i],rawy[i],det_id[i]]++
        endif
     endfor
     raw2det,fpd,fp
     writefits,obsid+ab[t]+'.'+key+'.raw.fits',fp

     table=readfits(files_opt[t],header,EXTEN_NO=1,/SILENT)
     time=TBGET( header, table, 'TIME', /NOSCALE)
     x_oa=TBGET( header, table, 'X_OA', /NOSCALE)
     y_oa=TBGET( header, table, 'Y_OA', /NOSCALE)

     for i=0L,N_ELEMENTS(x_oa)-1 do begin
        if(x_oa[i] ge 0 and y_oa[i] ge 0 and x_oa[i] lt skyw and y_oa[i] lt skyw) then begin
           optsky[x_oa[i],y_oa[i],t]++
        endif
     endfor
     
     ;;
     ;; Find a peak 
     ;;
     ii_max=0L
     jj_max=0L
     vv_max=-10L
     for ii=1L,skyw-2 do begin
        for jj=1L,skyw-2 do begin
           if(optsky[ii,jj,t] lt 100) then continue
           p00=optsky[ii-1,jj-1,t]
           p10=optsky[ii,jj-1,t]
           p20=optsky[ii+1,jj-1,t]

           p01=optsky[ii-1,jj,t]
           p11=optsky[ii,jj,t]
           p21=optsky[ii+1,jj,t]

           p02=optsky[ii-1,jj+1,t]
           p12=optsky[ii,jj+1,t]
           p22=optsky[ii+1,jj+1,t]

           if(p11 gt p00 and $
              p11 gt p10 and $
              p11 gt p20 and $
              p11 gt p01 and $
              p11 gt p21 and $
              p11 gt p02 and $
              p11 gt p12 and $
              p11 gt p12 and $
              p11 gt p22) then begin
              if(vv_max lt p11) then begin
                 vv_max=p11
                 ii_max=ii
                 jj_max=jj
              endif
           endif
        endfor
     endfor


     writefits,obsid+ab[t]+'.'+key+'.sky_oa.fits',optsky[*,*,t],hdr

     skydetfile=obsid+ab[t]+'.'+key+'.skydet.fits'
     xyad,hdr,ii_max,jj_max,ra,dec
     print,'OA peak position (RA/Dec): ',ra,dec
     print, sphdist(ra, dec, ra_obj, dec_obj, /DEGREES)*3600
     cmd='nuskytodet pntra='+String(RA_OBJ,Format='(f9.4)')+$
         ' pntdec='+String(DEC_OBJ,Format='(f9.4)')+$
         ' mastaspectfile='+files_mst+$
         ' attfile='+files_att+$
         ' clobber=yes chatter=0 skydetfile='+skydetfile+' instrument=FPM'+ab[t]
     spawn,cmd
     if(~file_test(skydetfile)) then begin
        print,'Does nuskytodet work?'
        stop
     endif
     table=readfits(skydetfile,header,EXTEN_NO=1,/SILENT)
     det1x=TBGET( header, table, 'DET1X', /NOSCALE)
     det1y=TBGET( header, table, 'DET1Y', /NOSCALE)
     nuplan_instrmap,ab[t],pixmapx=pixmapx,pixmapy=pixmapy,detnum=detnum
     fpd = lonarr(npix/2.,npix/2.,4)
     for i=0L,N_ELEMENTS(det1x)-1 do begin
        if(det1x[i] ge 0 and det1y[i] ge 0 and det1x[i] lt det1w and det1y[i] lt det1w) then begin
           optdet[det1x[i],det1y[i],t]++
           xpix=pixmapx[det1x[i],det1y[i]]
           ypix=pixmapy[det1x[i],det1y[i]]
           ndet=detnum[det1x[i],det1y[i]]
           if(xpix eq -1 or ypix eq -1) then continue
           fpd[xpix, ypix, ndet]++
        endif
     endfor

     raw2det,fpd,fp
     
     writefits,obsid+ab[t]+'.'+key+'.raw_oa.fits',fp
     writefits,obsid+ab[t]+'.'+key+'.det_oa.fits',optdet[*,*,t]

     ; see http://www.stsci.edu/~mperrin/software/sources/fwcentroid.pro
     sum=0.0
     XSUM = 0.0
     XSUM2 = 0.0
     XSUM3 = 0.0
     YSUM = 0.0
     YSUM2 = 0.0
     YSUM3 = 0.0
     for ii=0L,npix-2 do begin
        for jj=1L,npix-2 do begin
           if(fp[ii,jj] lt 100) then begin
              XSUM += ii * fp[ii,jj]
              XSUM2 += ii^2 * fp[ii,jj]
              XSUM3 += ii^3 * fp[ii,jj]
              YSUM += jj * fp[ii,jj]
              YSUM2 += jj^2 * fp[ii,jj]
              YSUM3 += jj^3 * fp[ii,jj]
              SUM += fp[ii,jj]
           endif
        endfor
     endfor
     XCEN = XSUM / SUM
     XMOMENT2 = XSUM2 / SUM
     XMOMENT3 = XSUM3 / SUM
     YCEN = YSUM / SUM
     YMOMENT2 = YSUM2 / SUM
     YMOMENT3 = YSUM3 / SUM

     print,ab[t],XCEN,YCEN,det2mm(XCEN),det2mm(YCEN)
     ;print,'2,-2: ',mm2det(2),mm2det(-2)
     

     oax_mm[t]=det2mm(XCEN)
     oay_mm[t]=det2mm(YCEN)
  endfor
  

end
