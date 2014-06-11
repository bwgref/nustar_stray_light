
; e.g.
;root='/Data/NuSTAR/data/mini-survey'
;obsid_array=['40010001002', '40010002001', '40010003001', '40010004001', '40010005001', '40010006001']

pro make_hotpix_raw,root=root,obsid=obsid_array,silent=silent,n_scales=n_scales, max_scale=max_scale, threshold=threshold

  if(n_elements(n_scales) eq 0) then n_scales=7
  if(n_scales le 0) then begin
     print,'Wrong n_scales: ',n_scales 
     stop
  endif
  if(n_elements(max_scale) eq 0) then max_scale=3
  if(max_scale le 0 or max_scale gt n_scales) then begin
     print,'Wrong max_scale: ',max_scale 
     stop
  endif
  if(n_elements(silent) eq 0) then silent=0
  if(n_elements(threshold) eq 0) then threshold=5.
  if(n_elements(border) eq 0) then border=20
  if(n_elements(root) eq 0) then begin
     print,"Give me root directory! root='path to data'" 
     return
  endif
  if(n_elements(obsid_array) eq 0) then begin
     print,"Give me obsid array! e.g. obsid=['40010001002', '40010002001',..]"
     return
  endif

  npix=64
  raw=lonarr(npix/2., npix/2., 4, 2)
  for io=0,n_elements(obsid_array)-1 do begin
     obsid=obsid_array(io)
        files_evt = root+'/'+obsid+['/event_cl/nu'+obsid+'A01_cl.evt','/event_cl/nu'+obsid+'B01_cl.evt']
        for t=0,1 do begin
           print,'Read ',files_evt[t]
           table=readfits(files_evt[t],header,EXTEN_NO=1,/SILENT)
           x=TBGET( header, table, 'RAWX', /NOSCALE)
           y=TBGET( header, table, 'RAWY', /NOSCALE)
           d=TBGET( header, table, 'DET_ID', /NOSCALE)
           for i=0L,N_ELEMENTS(x)-1 do begin
              if(x[i] gt 0 and y[i] gt 0) then raw[x[i],y[i],d[i],t]++
           endfor
        endfor
     endfor

  mod_name=['A','B']
  badpix_cnt=LONARR(2)
  hotpix=LONARR(npix,npix,2)
  det=LONARR(npix,npix,2)
  raw2det,raw[*,*,*,0],fpA
  raw2det,raw[*,*,*,1],fpB
  if(silent ne 1) then begin 
     writefits,'raw_stack_ima_fpA.img',fpA
     writefits,'raw_stack_ima_fpB.img',fpB
  endif
  det[*,*,0]=fpA
  det[*,*,1]=fpB
  
  for t=0,1 do begin
     map=LONARR(npix,npix)
     ATROUS, det[*,*,t], decomposition = dwt, n_scales=n_scales
     bkg=dwt[*,*,0]  
     for s=1,max_scale-1 do bkg+=dwt[*,*,s]     
     for i=0,npix-1 do begin
        for j=0,npix-1 do begin
           if(det[i,j,t] gt bkg[i,j]*threshold and $
              bkg[i,j] gt 0. and $
              i ge border and j ge border and $
              i le (npix-border-1) and j le (npix-border-1)) then begin
              map[i,j]=det[i,j,t]/bkg[i,j]
              hotpix[i,j]=1
              badpix_cnt[t]++
           endif
        endfor
     endfor

     if(silent ne 1) then begin 
        print,'fp'+mod_name[t]+': flag ',badpix_cnt[t]
        writefits,'raw_stack_ima_fp'+mod_name[t]+'.img',det[*,*,t]
        writefits,'raw_stack_bkg_fp'+mod_name[t]+'.img',bkg
        writefits,'raw_stack_excess_fp'+mod_name[t]+'.img',map
     endif

  endfor
  save,filename='raw_stack_hotpix.sav', hotpix, root, obsid_array

end
