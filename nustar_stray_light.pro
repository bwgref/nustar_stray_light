pro get_oa,oa,pos=pos
  common nuplan, nu, status, sources, target

   if(nu.oa[0] gt max(nu.xpos_array) or $
      nu.oa[0] lt min(nu.xpos_array) or $
      nu.oa[1] gt max(nu.ypos_array) or $
      nu.oa[1] lt min(nu.ypos_array)) then begin
      print,'Error: OA is out of FOV'
      widget_control, status.mainid, /destroy
      stop
   endif

  xpos_near = Min(Abs(nu.xpos_array - oa[0]), xpos_index)
  ypos_near = Min(Abs(nu.ypos_array - oa[1]), ypos_index)

  dx=1
  dy=1

  if(nu.xpos_array[xpos_index] gt oa[0]) then dx=-1
  if(nu.ypos_array[ypos_index] gt oa[1]) then dy=-1

  pos=fltarr(2)
  pos[0]=(oa[0]-nu.xpos_array[xpos_index+dx])/(nu.xpos_array[xpos_index]-nu.xpos_array[xpos_index+dx])*(-dx)+(xpos_index+dx)
  pos[1]=(oa[1]-nu.ypos_array[ypos_index+dy])/(nu.ypos_array[ypos_index]-nu.ypos_array[ypos_index+dy])*(-dy)+(ypos_index+dy)

end

pro plot_targets
  common nuplan, nu, status, sources, target

  n_detx=nu.n_detx
  n_dety=nu.n_dety
  dr=nu.dr
  rd=nu.rd

  make_astr_nustar, status.ra, status.dec, status.pa, astr=astr_cur, oa=mm2det(nu.oa)
  ;print,'shift: ',nu.fov_shift_x,nu.fov_shift_y
  xy2ad,mm2det(nu.oa[0]+nu.fov_shift_x)-1,mm2det(nu.oa[1]+nu.fov_shift_y)-1,astr_cur,new_ra,new_dec
  ;print,'New RA/Dec: ',new_ra,new_dec,status.ra, status.dec
  make_astr_nustar, new_ra, new_dec, status.pa, astr=astr, oa=mm2det(nu.oa)

   ;ad2xy,status.ra, status.dec, astr_cur, oax,oay
   ;print,'mm2det:',mm2det(nu.oa[0]),mm2det(nu.oa[1])
   ;print,'OA position and index',nu.oa[0],nu.oa[1],' ->',oax,oay

  ;status.ra=new_ra
  ;status.dec=new_dec

   green = GETCOLOR('green', 100)
   wt = GETCOLOR('white', 100)
   PLOTS, CIRCLE(nu.fov_shift_x, nu.fov_shift_y, 0.5), color=wt
   for jj=0, n_elements(target.src_ra)-1 do begin 
      if(target.src_flag[jj] eq 0) then continue
      ad2xy,target.src_ra[jj],target.src_dec[jj],astr,dx,dy
      if(dx gt 0 and dx lt n_detx and dy gt 0 and dy lt n_dety) then begin
         ;print,'TARGET:',target.src_name[jj],dx,dy,target.src_ra[jj],target.src_dec[jj],det2mm(dx),det2mm(dy)
         d=sphdist(target.src_ra[jj],target.src_dec[jj],new_ra,new_dec,/DEGREES)*60.
         print,'TARGET OA offset, arcmin: ',target.src_name[jj], d,format='(a,1x,a,1x,f8.2)'
      endif else begin
         continue
      endelse
      push,ptx,det2mm(dx+1)
      push,pty,det2mm(dy+1)
      push,ptn,target.src_name[jj]
   endfor
   
   for n=0,n_elements(ptx)-1 do begin
      PLOTS, CIRCLE(ptx[n], pty[n], 2.5), color=green
      ;oplot, [ptx[n]], [pty[n]], psym=1, color=green
      ;xyouts, ptx[n] , pty[n]+3.0, ptn[n], CHARSIZE=0.9, color=green
   endfor
  
end




pro nustar_stray_light, pnt_ra, pnt_dec, pa=pnt_pa , stray_catalog=stray_catalog, oa=oa, oaa=oaa, oab=oab, smooth=smooth, $
            silent=silent, quit=quit, target_catalog=target_catalog, scan_step=scan_step, badpix=badpix, $
            key=key,slp_level=slp_level, fmin=fmin, do_scan = do_scan
common nuplan, nu, status, sources, target
  forward_function read_combined_catalog

  if(n_elements(smooth) eq 0) then smooth=0  
  if(n_elements(key) eq 0) then key='nuplan'  
  if(n_elements(badpix) eq 0) then badpix=0  
  if(n_elements(oa) eq 0) then oa=[3., 3.] ; "default" position of OA in mm, will be converted to crpix WCS keyword
  if(n_elements(oaa) eq 0) then oaa=[3., 3.] 
  if(n_elements(oab) eq 0) then oab=[3., 3.] 
  if(n_elements(scan_step) eq 0) then scan_step=5.
  if(n_elements(silent) eq 0) then silent=0
  if(n_elements(quit) eq 0) then quit=0
  if n_elements(fmin) eq 0 then fmin = 5. ; default to 5 mCrab 
  if keyword_set(do_scan) then pnt_pa = 0

  target={src_name:[''],src_ra:[0.0],src_dec:[0.0],src_flag:[0],index:[1]}

  err = 0
  if(n_elements(slp_level) eq 0) then slp_level=0      
  if(n_elements(pnt_ra) eq 0) then err = 1 ;read,'NuPLAN: Enter R.A.: ',pnt_ra      
  if(n_elements(pnt_dec) eq 0) then err = 1;read,'NuPLAN: Enter Dec.: ',pnt_dec      
  if(n_elements(pnt_pa) eq 0) then err = 1;read,'NuPLAN: Enter PA: ',pnt_pa 
  ;if(n_elements(stray_catalog) eq 0) then read,'NuPLAN: Enter path to the source input catalog: ',stray_catalog 
  if err eq 1 then message, "Syntax is: nustar_stray_light, ra, dec, pa=PA or nustar_stray_light, ra, dec, /do_scan"


  ; Check for existing output files and remove
  outfiles = ['pa_scan.dat', 'scan.pdf', 'scan.ps']
  fi = file_info(outfiles)
  for i = 0, n_elements(outfiles) - 1 do begin
     if fi[i].exists then spawn, 'rm '+outfiles[i]
     ; Set up empty files in case things crash out below.
     spawn, 'touch '+outfiles[i]
  endfor


  hgap = 0.15                   ; detector half gap [mm]                                                                                                          
  n_detx = 64                   ; number of detector x bins 
  n_dety = 64                   ; number of detector y bins 
  d0mask = fltarr(n_detx, n_dety)
  d1mask = fltarr(n_detx, n_dety)
  xpos_array = fltarr(n_detx)
  ypos_array = fltarr(n_dety)
  xpos_2d_array=fltarr(n_detx,n_dety)
  ypos_2d_array=fltarr(n_detx,n_dety)
  for xx = 0, n_detx-1 do begin
     for yy = 0, n_dety-1 do begin
        q_x = sign(xx+1 - 32.5)                                        
        q_y = sign(yy+1 - 32.5)
        xpos_array(xx) = det2mm(xx)
        ypos_array(yy) = det2mm(yy)
      ;xpos_array(xx) = (xx+1-32.5)*0.6048 + hgap*sign(xx+1-32.5)
      ;ypos_array(yy) = (yy+1-32.5)*0.6048 + hgap*sign(yy+1-32.5)
        xpos_2d_array(xx,yy)=xpos_array(xx)
        ypos_2d_array(xx,yy)=ypos_array(yy)
     endfor
  endfor


  status={ra:pnt_ra, dec:pnt_dec, ra_orig:pnt_ra, dec_orig:pnt_dec, pa:pnt_pa,$
          slpa:0,slpb:0,eff:0,vis:0,winid:0,mainid:0,infoid:0,TopTab:0,table:0,$
          infotext:0,silent:silent,scan_step:scan_step, key:key, smooth:smooth, $
          loss0:0., loss1:0., ghost_ray:0.}
  nu={hgap:hgap,n_detx:n_detx,n_dety:n_dety,xpos_array:xpos_array,ypos_array:ypos_array,oa:oa,oaa:oaa,oab:oab,oa_prev:oa,fov_shift_x:0.0D,fov_shift_y:0.0D,fov_shift_step:0.5D,dr:!PI/180.,rd:180./!PI}

; Read in the combined BAT 70 Month and INTEGRAL galactic plane survey catalogs
;  read_combined_catalog, ra=pnt_ra, dec=pnt_dec, src_name, src_ra, src_dec, src_flux, src_flag, $
;                         fmin=fmin
  read_3to30_catalog, ra = pnt_ra, dec =pnt_dec, src_name, src_ra, src_dec, src_flux, src_flag, $
                      int_src, fmin = fmin
  sources={src_name:src_name,src_ra:src_ra,src_dec:src_dec,src_flux:src_flux,src_flag:src_flag,int_src:int_src}


  if ~keyword_set(do_scan) then begin
     stray_light_render, badpix=badpix
  endif else begin

;   print,'Running PA scan 0,360,5 deg'
     openw, lun,  'pa_scan.dat', /get_lun
                                ; save source list used in this run
;     for ii=0, n_elements(sources.src_name)-1 do printf,lun,'# ',$
;        sources.src_name(ii),',',sources.src_ra(ii),', ',sources.src_dec(ii),',',sources.src_flag(ii)
     save_pa=status.pa
     save_silent=status.silent
     status.silent=1

     setps
     printf, lun, 'PA FPMA_LOSS FPMB_LOSS FPMA_DET0_LOSS FPMB_DET0_LOSS'
     for s=0.0,360.0,status.scan_step do begin
        status.pa=s
        stray_light_render
;      print, status.pa, status.slpa, status.slpb, status.loss0, status.loss1,format='(5f8.2)'
        printf,lun, status.pa, status.slpa, status.slpb, status.loss0, status.loss1,format='(5f8.2)'
     
     endfor
     close, lun
     free_lun, lun
     endps
     spawn, 'mv idlout.pdf scan.pdf'
     status.pa=save_pa
     status.silent=save_silent
  endelse


end
