pro map_hotpix,root=root,obsid=obsid_array,silent=silent

  if(n_elements(silent) eq 0) then silent=0
  if(n_elements(threshold) eq 0) then threshold=5.
  if(n_elements(root) eq 0) then begin
     print,'Give me root directory! root=\'dir to data\' exit.' 
     return
  endif
  if(n_elements(obsid_array) eq 0) then begin
     print,'Give me obsid array! obsid=[\'40010001002\', \'40010002001\',..] exit.' 
     return
  endif

  det1w=360
  border=20

  ;root='/Data/NuSTAR/data/mini-survey'
  ;obsid_array=['40010001002', '40010002001', '40010003001', '40010004001', '40010005001', '40010006001']


  det1=LONARR(det1w,det1w,2)
  for io=0,n_elements(obsid_array)-1 do begin
     obsid=obsid_array(io)
        files_evt = root+'/'+obsid+['/event_cl/nu'+obsid+'A01_cl.evt','/event_cl/nu'+obsid+'B01_cl.evt']
        for t=0,1 do begin
           print,'Read ',files_evt[t]
           table=readfits(files_evt[t],header,EXTEN_NO=1,/SILENT)
           x=TBGET( header, table, 'DET1X', /NOSCALE)
           y=TBGET( header, table, 'DET1Y', /NOSCALE)
           for i=0L,N_ELEMENTS(x)-1 do begin
              if(x[i] gt 0 and y[i] gt 0) then det1[x[i],y[i],t]++
           endfor
        endfor
     endfor

  mod_name=['A','B']
  badpix_cnt=LONARR(2)
  for t=0,1 do begin
     map=LONARR(det1w,det1w)
     ATROUS, det1[*,*,t], decomposition = dwt, n_scales=6
     bkg=dwt[*,*,0]+dwt[*,*,1]+dwt[*,*,2];+dwt[*,*,3];+dwt[*,*,4]
     ;bkg_smo = filter_image( bkg, SMOOTH=10)
     
     for i=0,det1w-1 do begin
        for j=0,det1w-1 do begin
           if(det1[i,j,t] gt bkg[i,j]*threshold and $
              bkg[i,j] gt 0. and $
              i ge border and j ge border and $
              i le (det1w-border-1) and j le (det1w-border-1)) then begin
              map[i,j]=det1[i,j,t]/bkg[i,j]
              badpix_cnt[t]++
           endif
        endfor
     endfor
     print,'fp'+mod_name[t]+': flag ',badpix_cnt[t]
     writefits,'det1_stack_ima_fp'+mod_name[t]+'.img',det1[*,*,t]
     writefits,'det1_stack_bkg_fp'+mod_name[t]+'.img',bkg
     writefits,'det1_stack_excess_fp'+mod_name[t]+'.img',map
  endfor

end
