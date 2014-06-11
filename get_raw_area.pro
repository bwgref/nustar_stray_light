
; e.g.
;root='/Data/NuSTAR/data/mini-survey'
;obsid_array=['40010001002', '40010002001', '40010003001', '40010004001', '40010005001', '40010006001']

pro get_raw_area, root=root, obsid=obsid, silent=silent, eventsdir=eventsdir, area=area, mask_det1=mask_det1, mask_raw=mask_raw, mask_sky=mask_sky
  det1w=360
  hotpix_max=1000
  ;e1=2.
  ;e2=10.
  ;cha1=(e1-1.6)/0.04
  ;cha2=(e2-1.6)/0.04
  if(n_elements(silent) eq 0) then silent=0
  if(n_elements(eventsdir) eq 0) then eventsdir='event_cl'
  if(n_elements(root) eq 0) then begin
     print,"Give me root directory! root='path to data'" 
     return
  endif
  if(n_elements(obsid) eq 0) then begin
     print,"Give me obsid! e.g. obsid='40010001002'"
     return
  endif

  npix=64
  mask_raw=lonarr(npix,npix,2)  
  badpix_raw=lonarr(npix,npix,2)  
  fpdA = lonarr(npix/2.,npix/2.,4) 
  fpdB = lonarr(npix/2.,npix/2.,4) 

  mask_fpdA = lonarr(npix/2.,npix/2.,4)+1 
  mask_fpdB = lonarr(npix/2.,npix/2.,4)+1

  mod_name=['A','B']
  mask_det1=LONARR(det1w,det1w,2)
  files_evt = root+'/'+obsid+'/'+eventsdir+['/nu'+obsid+'A01_cl.evt','/nu'+obsid+'B01_cl.evt']
  for t=0,1 do begin
     print,'Read ',files_evt[t]

     table=readfits(files_evt[t],header,EXTEN_NO=1,/SILENT)
     ;grade=TBGET( header, table, 'GRADE', /NOSCALE)
     ;chan=TBGET( header, table, 'PI', /NOSCALE)
     ;det_id=TBGET( header, table, 'DET_ID', /NOSCALE)
     ;det1x=TBGET( header, table, 'DET1X', /NOSCALE)-1
     ;det1y=TBGET( header, table, 'DET1Y', /NOSCALE)-1
     rawx=TBGET( header, table, 'RAWX', /NOSCALE)
     rawy=TBGET( header, table, 'RAWY', /NOSCALE)

     ;for i=0L,N_ELEMENTS(det1x)-1 do begin
     ;   if(det1x[i] gt 0 and det1y[i] gt 0) then begin
     ;      det1[det1x[i],det1y[i],t]++
     ;   endif
     ;endfor
      
     for id=0,3 do begin
        table=readfits(files_evt[t],header,EXTEN_NO=(3+id),/SILENT)
        badrawx=TBGET( header, table, 'RAWX', /NOSCALE)
        badrawy=TBGET( header, table, 'RAWY', /NOSCALE)
        for r=0,n_elements(badrawx)-1 do begin
           bx=badrawx[r]
           by=badrawy[r]
           if(t eq 0) then begin
              fpdA[bx,by,id]++
              mask_fpdA[bx,by,id]=0
           endif else begin
              fpdB[bx,by,id]++
              mask_fpdB[bx,by,id]=0
           endelse
        endfor
     endfor
  endfor

  mask_fpdA[0,*,*]=0
  mask_fpdA[31,*,*]=0
  mask_fpdA[*,0,*]=0
  mask_fpdA[*,31,*]=0

  mask_fpdB[0,*,*]=0
  mask_fpdB[31,*,*]=0
  mask_fpdB[*,0,*]=0
  mask_fpdB[*,31,*]=0

  fpdA[0,*,*]=1
  fpdA[31,*,*]=1
  fpdA[*,0,*]=1
  fpdA[*,31,*]=1

  fpdB[0,*,*]=1
  fpdB[31,*,*]=1
  fpdB[*,0,*]=1
  fpdB[*,31,*]=1

  raw2det,fpdA,fpA
  raw2det,fpdB,fpB

  raw2det,mask_fpdA,mask_fpA
  raw2det,mask_fpdB,mask_fpB
  mask_raw[*,*,0]=mask_fpA
  mask_raw[*,*,1]=mask_fpB
  
  ;if(silent ne 1) then begin 
  ;   writefits,'raw_badpix_fpA.img',fpA
  ;   writefits,'raw_badpix_fpB.img',fpB
  ;endif

  badpix_raw[*,*,0]=fpA
  badpix_raw[*,*,1]=fpB

  totA=n_elements(where(fpA eq 0))
  totB=n_elements(where(fpB eq 0))
  area=[totA,totB]*(0.6*0.1)^2 ; 0.6mm
  
  ; RAW -> DET1
  mask_det1[*]=0
  for t=0,1 do begin
     nuplan_instrmap,mod_name[t],pixmapx=pixmapx,pixmapy=pixmapy,detnum=detnum
     for i=0L,det1w-1 do for j=0L,det1w-1 do begin
        if(pixmapx[i,j] eq -1 or pixmapy[i,j] eq -1) then continue
        if(t eq 0) then begin
           if (fpdA[pixmapx[i,j], pixmapy[i,j],detnum[i,j]] eq 0) then mask_det1[i,j,t]=1
        endif
        if(t eq 1) then begin
           if (fpdB[pixmapx[i,j], pixmapy[i,j],detnum[i,j]] eq 0) then mask_det1[i,j,t]=1
        endif
     endfor
  endfor

     ;writefits,'det1_mask_fpA.img', mask_det1[*,*,0]
     ;writefits,'det1_mask_fpB.img', mask_det1[*,*,1]
     ;writefits,'raw_mask_fpA.img', mask_raw[*,*,0]
     ;writefits,'raw_mask_fpB.img', mask_raw[*,*,1]
     

end
