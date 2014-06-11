
; e.g.
;root='/Data/NuSTAR/data/mini-survey'
;obsid_array=['40010001002', '40010002001', '40010003001', '40010004001', '40010005001', '40010006001']

pro make_det1,root=root,obsid=obsid_array,postfix=postfix,silent=silent, grade_max=grade_max, e1=e1, e2=e2, key=key, livetime=livetime, evtpath=evtpath, det1bad=det1bad_map, mode=mode

  npix=64
  det1w=360
  half=180
  gapx1=[175,177]
  gapx2=[189,191]
  gapy1=[175,177]
  gapy2=[189,191]
  hotpix_max=10000
  ab=['A','B']

  if(n_elements(e1) eq 0) then e1=3.
  if(n_elements(e2) eq 0) then e2=10.
  cha1=(e1-1.6)/0.04
  cha2=(e2-1.6)/0.04

  if(n_elements(evtpath) eq 0) then evtpath='event_cl'
  if(n_elements(key) eq 0) then key='key'
  if(n_elements(postfix) eq 0) then postfix=''
  if(n_elements(mode) eq 0) then mode='01'
  if(n_elements(grade_max) eq 0) then grade_max=26
  if(n_elements(silent) eq 0) then silent=0
  if(n_elements(root) eq 0) then begin
     print,"Give me root directory! root='path to data'" 
     return
  endif
  if(n_elements(obsid_array) eq 0) then begin
     print,"Give me obsid array! e.g. obsid=['40010001002', '40010002001',..]"
     return
  endif

  stack_raw=lonarr(npix,npix,2)
  stack_det1=LONARR(det1w,det1w,2)
  livetime_total=DBLARR(2)
  
  chip_rate=LONARR(4,2)
  for io=0,n_elements(obsid_array)-1 do begin
     obsid=obsid_array(io)
     stack_fpdA = lonarr(npix/2.,npix/2.,4) 
     stack_fpdB = lonarr(npix/2.,npix/2.,4) 
     det1=LONARR(det1w,det1w,2)
     det1bad_map=LONARR(det1w,det1w,2)+1
     livetime=DBLARR(2)
     
     

     files_evt = root+'/'+obsid+['/'+evtpath+'/nu'+obsid+'A'+mode+'_cl.evt','/'+evtpath+'/nu'+obsid+'B'+mode+'_cl.evt']+postfix
     get_raw_area,root=root,obsid=obsid,/silent,eventsdir=evtpath,area=area,mask_det1=mask_det1,mask_raw=mask_raw
     for t=0,1 do begin
        mkhdr, hdr_long, 3, [det1w,det1w]
        SXADDPAR, hdr_long, 'AREA', area[t], ' Area, sq. cm'
        SXADDPAR, hdr_long, 'OBSID', obsid, ''
        if(silent ne 1) then begin
           writefits,obsid+ab[t]+'.'+key+'.raw_msk.fits',mask_raw[*,*,t],hdr_long
           writefits,obsid+ab[t]+'.'+key+'.det_msk.fits',mask_det1[*,*,t],hdr_long
        endif
     endfor
     
     
     for t=0,1 do begin
        if(~file_test(files_evt[t])) then continue
        print,'Read ',files_evt[t]
        table=readfits(files_evt[t],header,EXTEN_NO=1,/SILENT)
        livetime[t] = SXPAR(header,'LIVETIME',COUNT=A)
        grade=TBGET( header, table, 'GRADE', /NOSCALE)
        chan=TBGET( header, table, 'PI', /NOSCALE)
        det_id=TBGET( header, table, 'DET_ID', /NOSCALE)
        det1x=TBGET( header, table, 'DET1X', /NOSCALE)-1
        det1y=TBGET( header, table, 'DET1Y', /NOSCALE)-1
        rawx=TBGET( header, table, 'RAWX', /NOSCALE)
        rawy=TBGET( header, table, 'RAWY', /NOSCALE)
        index=where(grade ge 0 and grade le grade_max and chan gt cha1 and chan lt cha2, count)
        if (count eq 0) then continue
        print, 'Read events: ', count, n_elements(grade), ' livetime: ',livetime[t]
        grade=grade(index)
        chan=chan(index)
        det_id=det_id(index)
        det1x=det1x(index)
        det1y=det1y(index)
        rawx=rawx(index)
        rawy=rawy(index)

	livetime_total[t]+=livetime[t]

        mkhdr, hdr_det1_long, 3, [det1w,det1w]
        SXADDPAR, hdr_det1_long, 'AREA', area[t], ' Area, sq. cm'
        SXADDPAR, hdr_det1_long, 'LIVETIME', livetime[t], ' Exposure, s.'
        SXADDPAR, hdr_det1_long, 'EMIN', e1, ' EMIN, keV'
        SXADDPAR, hdr_det1_long, 'EMAX', e2, ' EMIN, keV'
        SXADDPAR, hdr_det1_long, 'OBSID', obsid, ''
        SXADDPAR, hdr_det1_long, 'MODULE', ab[t], ' Module'

        mkhdr, hdr_raw_long, 3, [npix,npix]
        SXADDPAR, hdr_raw_long, 'AREA', area[t], ' Area, sq. cm'
        SXADDPAR, hdr_raw_long, 'LIVETIME', livetime[t], ' Exposure, s.'
        SXADDPAR, hdr_raw_long, 'EMIN', e1, ' EMIN, keV'
        SXADDPAR, hdr_raw_long, 'EMAX', e2, ' EMIN, keV'
        SXADDPAR, hdr_raw_long, 'OBSID', obsid, ' Obs. ID'
        SXADDPAR, hdr_raw_long, 'MODULE', ab[t], ' Module'
        
        for i=0L,N_ELEMENTS(det_id)-1 do begin
           if(mask_raw[rawx[i],rawy[i],t] eq 0) then continue
           chip_rate[det_id[i],t]++
           if(t eq 0) then stack_fpdA[rawx[i],rawy[i],det_id[i]]++
           if(t eq 1) then stack_fpdB[rawx[i],rawy[i],det_id[i]]++
        endfor
        
        for i=0L,N_ELEMENTS(det1x)-1 do begin
           if(det1x[i] ge 0 and det1y[i] ge 0) then begin
              if(mask_det1[det1x[i],det1y[i],t] eq 1) then begin 
                 det1[det1x[i],det1y[i],t]++
                 stack_det1[det1x[i],det1y[i],t]++
              endif
           endif
        endfor
        
        raw2det,stack_fpdA,stack_fpA
        raw2det,stack_fpdB,stack_fpB
        if(silent ne 1) then begin
           writefits,obsid+ab[t]+'.'+key+'.raw_ima.fits',(t eq 0)? stack_fpA: stack_fpB,hdr_det1_long
           writefits,obsid+ab[t]+'.'+key+'.det_ima.fits',det1[*,*,t],hdr_det1_long
        endif
     endfor
     endfor

  print,'Total A',total(stack_det1[*,*,0])
  print,'Total B',total(stack_det1[*,*,1])

        mkhdr, hdr_det1, 3, [det1w,det1w]
        SXADDPAR, hdr_det1, 'LIVETIME', livetime_total[0], ' Exposure, s.'
        SXADDPAR, hdr_det1, 'EMIN', e1, ' EMIN, keV'
        SXADDPAR, hdr_det1, 'EMAX', e2, ' EMIN, keV'
  writefits,'stack'+mode+'_fpA.'+key+'.det1.fits',stack_det1[*,*,0],hdr_det1_long
        SXADDPAR, hdr_det1, 'LIVETIME', livetime_total[1], ' Exposure, s.'
  writefits,'stack'+mode+'_fpB.'+key+'.det1.fits',stack_det1[*,*,1],hdr_det1_long



end
