
; e.g.
;root='/Data/NuSTAR/data/mini-survey'
;obsid_array=['40010001002', '40010002001', '40010003001', '40010004001', '40010005001', '40010006001']

pro make_hotpix_det1,root=root,obsid=obsid_array,silent=silent,n_scales=n_scales, max_scale=max_scale, grade_max=grade_max, e1=e1, e2=e2, key=key, livetime=livetime, evtpath=evtpath, det1bad=det1bad_map
  npix=64
  det1w=360
  half=180
  gapx1=[175,177]
  gapx2=[189,191]
  gapy1=[175,177]
  gapy2=[189,191]
  hotpix_max=10000

  if(n_elements(e1) eq 0) then e1=3.
  if(n_elements(e2) eq 0) then e2=10.
  cha1=(e1-1.6)/0.04
  cha2=(e2-1.6)/0.04

  if(n_elements(evtpath) eq 0) then evtpath='event_cl'
  if(n_elements(key) eq 0) then key='key'
  if(n_elements(n_scales) eq 0) then n_scales=7
  if(n_elements(grade_max) eq 0) then grade_max=4
  if(n_scales le 0 or n_scales gt 9) then begin
     print,'Wrong n_scales: ',n_scales 
     stop
  endif
  if(n_elements(max_scale) eq 0) then max_scale=3
  if(max_scale le 0 or max_scale gt n_scales) then begin
     print,'Wrong max_scale: ',max_scale 
     stop
  endif
  if(n_elements(silent) eq 0) then silent=0
  if(n_elements(threshold) eq 0) then threshold=10.
  if(n_elements(border) eq 0) then border=10
  if(n_elements(root) eq 0) then begin
     print,"Give me root directory! root='path to data'" 
     return
  endif
  if(n_elements(obsid_array) eq 0) then begin
     print,"Give me obsid array! e.g. obsid=['40010001002', '40010002001',..]"
     return
  endif

  ;stack_raw=lonarr(npix,npix,2)
  stack_fpdA = lonarr(npix/2.,npix/2.,4) 
  stack_fpdB = lonarr(npix/2.,npix/2.,4) 

  chip_rate=LONARR(4,2)
  det1=LONARR(det1w,det1w,2)
  det1_noedge=LONARR(det1w,det1w,2)
  det1bad_map=LONARR(det1w,det1w,2)+1
  expo=DBLARR(det1w,det1w,2)
  livetime=DBLARR(2)
  for io=0,n_elements(obsid_array)-1 do begin
     obsid=obsid_array(io)

     
     files_evt = root+'/'+obsid+['/'+evtpath+'/nu'+obsid+'A01_cl.evt','/'+evtpath+'/nu'+obsid+'B01_cl.evt']
     files_evt_orig = root+'/'+obsid+['/'+evtpath+'/nu'+obsid+'A01_cl_orig.evt','/'+evtpath+'/nu'+obsid+'B01_cl_orig.evt']
        for t=0,1 do begin
           if(file_test(files_evt[t]) and ~file_test(files_evt_orig[t])) then begin
              ;print,'copy ',files_evt[t],' -> ',files_evt_orig[t]
              file_copy,files_evt[t],files_evt_orig[t]
           endif
        endfor
        
        for t=0,1 do begin
           if(~file_test(files_evt_orig[t])) then continue
           print,'Read ',files_evt[t]
           table=readfits(files_evt_orig[t],header,EXTEN_NO=1,/SILENT)
           ltime = SXPAR(header,'LIVETIME',COUNT=A)
           livetime[t]+=ltime
           grade=TBGET( header, table, 'GRADE', /NOSCALE)
           chan=TBGET( header, table, 'PI', /NOSCALE)
           det_id=TBGET( header, table, 'DET_ID', /NOSCALE)
           det1x=TBGET( header, table, 'DET1X', /NOSCALE)-1
           det1y=TBGET( header, table, 'DET1Y', /NOSCALE)-1
           rawx=TBGET( header, table, 'RAWX', /NOSCALE)
           rawy=TBGET( header, table, 'RAWY', /NOSCALE)
           index=where(grade ge 0 and grade le grade_max and chan gt cha1 and chan lt cha2, count)
           if (count eq 0) then stop
           print, 'Read 0-4 grade events: ', count, n_elements(grade), ' livetime: ',ltime
           grade=grade(index)
           chan=chan(index)
           det_id=det_id(index)
           det1x=det1x(index)
           det1y=det1y(index)
           rawx=rawx(index)
           rawy=rawy(index)

           for i=0L,N_ELEMENTS(det_id)-1 do begin
              if(rawx[i] eq 0 or rawx[i] eq 31) then continue
              if(rawy[i] eq 0 or rawy[i] eq 31) then continue
              chip_rate[det_id[i],t]++
              if(t eq 0) then stack_fpdA[rawx[i],rawy[i],det_id[i]]++
              if(t eq 1) then stack_fpdB[rawx[i],rawy[i],det_id[i]]++
           endfor


           for i=0L,N_ELEMENTS(det1x)-1 do begin
              if(det1x[i] gt 0 and det1y[i] gt 0) then begin
                 det1[det1x[i],det1y[i],t]++
                 if(rawx[i] eq 0 or rawx[i] eq 31 or rawy[i] eq 0 or rawy[i] eq 31) then begin
                    det1bad_map[det1x[i],det1y[i],t]=0
                 endif else begin
                    det1_noedge[det1x[i],det1y[i],t]++
                 endelse
              endif
           endfor

           for ii=0,det1w-1 do begin
              for jj=0,det1w-1 do begin
                 if(det1(ii,jj,t) gt 0) then expo[ii,jj,t]+=ltime
              endfor
           endfor


           if(t eq 0) then begin
              push, det1xA,det1x
              push, det1yA,det1y
              push, det_idA,det_id
              push, rawxA,rawx
              push, rawyA,rawy
           endif 

           if(t eq 1) then begin
              push, det1xB,det1x
              push, det1yB,det1y
              push, det_idB,det_id
              push, rawxB,rawx
              push, rawyB,rawy
           endif 

           get_raw_area,root=root,obsid=obsid,/silent,eventsdir=evtpath,area=area,mask_det1=det1bad,mask_raw=mask_raw
           mask_fpA=mask_raw[*,0]
           mask_fpB=mask_raw[*,1]
           det1bad_map=det1bad_map*det1bad
     

        endfor
     endfor

  raw2det,stack_fpdA,stack_fpA
  raw2det,stack_fpdB,stack_fpB
  writefits,'stack_raw_fpA.'+key+'.img',stack_fpA
  writefits,'stack_raw_fpB.'+key+'.img',stack_fpB
  writefits,'stack_raw_mask_fpA.'+key+'.img',mask_fpA
  writefits,'stack_raw_mask_fpB.'+key+'.img',mask_fpB


  mod_name=['A','B']
  badpix_cnt=LONARR(2)
  hotpix_det1=LONARR(det1w,det1w,2)
  badpix_list=LONARR(hotpix_max,2,2)
  print
  for t=0,1 do begin
     tt=avg(chip_rate[*,t])
     print,mod_name[t]+' chip rates:',chip_rate[0,t],chip_rate[1,t],chip_rate[2,t],chip_rate[3,t]
     print,mod_name[t]+' chip coeff:',chip_rate[0,t]/tt,chip_rate[1,t]/tt,chip_rate[2,t]/tt,chip_rate[3,t]/tt
  endfor
  print
  for t=0,1 do begin
     map=LONARR(det1w,det1w)
     ATROUS, det1[*,*,t], decomposition = dwt, n_scales=n_scales
     bkg=dwt[*,*,0]  
     for s=1,max_scale-1 do bkg+=dwt[*,*,s]     
     for i=0,det1w-1 do begin
        for j=0,det1w-1 do begin
           if(det1[i,j,t] gt bkg[i,j]*threshold and $
              bkg[i,j] gt 0. and $
              i ge border and j ge border and $
              i le (det1w-border-1) and j le (det1w-border-1)) then begin
              map[i,j]=det1[i,j,t]/bkg[i,j]
              hotpix_det1[i,j,t]+=det1[i,j,t]
              badpix_list[badpix_cnt[t],0,t]=i
              badpix_list[badpix_cnt[t],1,t]=j
              badpix_cnt[t]++
              if(badpix_cnt[t] ge hotpix_max) then stop              
           endif
        endfor
     endfor


     if(silent ne 1) then begin 
        map_clean=LONARR(det1w,det1w)
        map_clean_noedge=LONARR(det1w,det1w)
        map_rate=DBLARR(det1w,det1w)
        for ii=0,det1w-1 do begin
           for gg=gapx1[t]-1,gapx2[t]-1 do begin
              det1bad_map[gg,ii,t]=0
           endfor
           for gg=gapy1[t]-1,gapy2[t]-1 do begin
              det1bad_map[ii,gg,t]=0
           endfor

           for jj=0,det1w-1 do begin
              if(map(ii,jj) gt 0) then continue
              map_clean[ii,jj]=det1[ii,jj,t]
              map_clean_noedge[ii,jj]=det1_noedge[ii,jj,t]
              if(expo[ii,jj,t] gt 0.0) then map_rate[ii,jj]=det1[ii,jj,t]/expo[ii,jj,t]
           endfor
        endfor
        mkhdr, hdr_long, 3, [det1w,det1w]
        SXADDPAR, hdr_long, 'LIVETIME', livetime[t], ' Exposure, s.'
        SXADDPAR, hdr_long, 'EMIN', e1, ' EMIN, keV'
        SXADDPAR, hdr_long, 'EMAX', e2, ' EMIN, keV'
        print,'fp'+mod_name[t]+': flag ',badpix_cnt[t]
        writefits,'stack_det1_bad_fp'+mod_name[t]+'.'+key+'.img',det1bad_map[*,*,t],hdr_long
        writefits,'stack_det1_ima_fp'+mod_name[t]+'.'+key+'.img',det1[*,*,t],hdr_long
        writefits,'stack_det1_bkg_fp'+mod_name[t]+'.'+key+'.img',bkg
        writefits,'stack_det1_hotpix_fp'+mod_name[t]+'.'+key+'.img',map
        writefits,'stack_det1_ima_cl_fp'+mod_name[t]+'.'+key+'.img',map_clean,hdr_long
        writefits,'stack_det1_ima_cl_noedge_fp'+mod_name[t]+'.'+key+'.img',map_clean_noedge,hdr_long
        writefits,'stack_det1_exp_cl_fp'+mod_name[t]+'.'+key+'.img',expo[*,*,1]
        writefits,'stack_det1_rate_cl_fp'+mod_name[t]+'.'+key+'.img',map_rate
        ;for s=0,n_scales-1 do begin
        ;   writefits,'stack_det1_atrous'+string(s,format='(i1)')+'_fp'+mod_name[t]+'.'+key+'.img',dwt[*,*,s]
        ;endfor
     endif
  endfor

  hotpix_raw=lonarr(npix,npix,2)
  
  fpdA = lonarr(npix/2.,npix/2.,4) 
  fpdB = lonarr(npix/2.,npix/2.,4) 

  for t=0,1 do begin
     for k=0, badpix_cnt[t]-1 do begin
        ; make A
        if(t eq 0) then begin
           index=where (det1xA eq badpix_list[k,0,t] and det1yA eq badpix_list[k,1,t], count)
           arx=lonarr(npix/2.)
           ary=lonarr(npix/2.)
           ard=lonarr(npix/2.)
           for i=0L,count-1 do begin
              arx[rawxA[index[i]]]++   
              ary[rawyA[index[i]]]++   
              ard[det_idA[index[i]]]++   
           endfor
           isx=reverse(sort(arx))
           isy=reverse(sort(ary))
           isd=reverse(sort(ard))
           ;if(silent ne 1) then print,'A DET1',badpix_list[k,0,t]+1,badpix_list[k,1,t]+1,' -> RAW', $
           ;                           isx[0],isy[0],isd[0],' (',count,')', $
           ;                           format='(a,2i4,a,3i3,a,i4,a)'
           fpdA[isx[0],isy[0],isd[0]]++
        endif

        ; make B
        if(t eq 1) then begin
           index=where (det1xB eq badpix_list[k,0,t] and det1yB eq badpix_list[k,1,t], count)
           arx=lonarr(npix/2.)
           ary=lonarr(npix/2.)
           ard=lonarr(npix/2.)
           for i=0L,count-1 do begin
              arx[rawxB[index[i]]]++   
              ary[rawyB[index[i]]]++   
              ard[det_idB[index[i]]]++   
           endfor
           isx=reverse(sort(arx))
           isy=reverse(sort(ary))
           isd=reverse(sort(ard))
           ;if(silent ne 1) then print,'B DET1',badpix_list[k,0,t]+1,badpix_list[k,1,t]+1,' -> RAW', $
           ;                           isx[0],isy[0],isd[0],' (',count,')', $
           ;                           format='(a,2i4,a,3i3,a,i4,a)'
           fpdB[isx[0],isy[0],isd[0]]++
        endif
     endfor     
  endfor    
  raw2det,fpdA,fpA
  raw2det,fpdB,fpB
  if(silent ne 1) then begin 
     writefits,'stack_raw_hotpix_fpA.'+key+'.img',fpA
     writefits,'stack_raw_hotpix_fpB.'+key+'.img',fpB
  endif
  hotpix_raw[*,*,0]=fpA
  hotpix_raw[*,*,1]=fpB

  save, filename='stack_hotpix.sav', hotpix_det1, hotpix_raw, root, obsid_array

  ;;
  ;; old style, trying converting DET1->RAW
  ;;
  ;; hotpix_raw=lonarr(npix,npix,2)
  ;; ima_raw=lonarr(npix,npix,2)
  ;; for t=0,1 do begin
  ;;    for i=0,det1w-1 do begin
  ;;       for j=0,det1w-1 do begin
  ;;          ii=round((i-22)/5.)
  ;;          jj=round((j-22)/5.)
  ;;          if(ii ge 1 and jj ge 1 and ii lt npix-1 and jj lt npix-1) then begin
  ;;             if(hotpix_det1[i,j,t] eq 1) then begin
  ;;                hotpix_raw[ii,jj,t]=1
  ;;             endif
  ;;             ima_raw[ii,jj,t]+=det1[i,j,t]
  ;;          endif
  ;;       endfor
  ;;    endfor
  ;;    if(silent ne 1) then begin 
  ;;       writefits,'stack_raw_ima_fp'+mod_name[t]+'.img',ima_raw[*,*,t]
  ;;       writefits,'stack_raw_hotpix_fp'+mod_name[t]+'.img',hotpix_raw[*,*,t]
  ;;    endif
  ;; endfor
  ;; save, filename='stack_hotpix.sav', hotpix_det1, hotpix_raw, root, obsid_array

end
