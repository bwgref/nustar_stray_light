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

pro read_table
  common nuplan, nu, status, sources, target
  WIDGET_CONTROL, status.table, GET_VALUE=tab
  nrows=n_elements(tab)
  for i=0,nrows-1 do begin
     sources.src_name[i]=tab[i].name
     sources.src_ra[i]=tab[i].ra
     sources.src_dec[i]=tab[i].dec
     sources.src_flag[i]=tab[i].flag
  end
end

PRO nuplan_event, ev                                          ; event handler
  common nuplan, nu, status, sources, target
                                ;this common gives us access to variables in the main program
  widget_control, ev.id, get_uvalue=uvalue ; get the uvalue

  read_table

  numTabs = WIDGET_INFO(status.TopTab, /TAB_NUMBER) 
  thisTab = WIDGET_INFO(status.TopTab, /TAB_CURRENT) 
 
   green = GETCOLOR('green', 100)
   wt = GETCOLOR('white', 100)
   blue = GETCOLOR('blue', 100)

  if(n_elements(uvalue) eq 0) then return
  CASE uvalue OF                           ; choose case
     'pa' : begin
        status.pa=ev.value
     end
     'scan' : begin
        print,'Running PA scan 0,360,5 deg'
        openw, lun,  'pa_scan.dat', /get_lun
        ; save source list used in this run
        for ii=0, n_elements(sources.src_name)-1 do printf,lun,'# ',$
           sources.src_name(ii),',',sources.src_ra(ii),', ',sources.src_dec(ii),',',sources.src_flag(ii)
        save_pa=status.pa
        save_silent=status.silent
        status.silent=1
        for s=0.0,360.0,status.scan_step do begin
           status.pa=s
           stray_light_render
           print, status.pa, status.slpa, status.slpb, status.eff, status.vis,format='(5f8.2)'
           printf,lun, status.pa, status.slpa, status.slpb, status.eff, status.vis,format='(5f8.2)'
        endfor
        close, lun
        free_lun, lun
        status.pa=save_pa
        status.silent=save_silent
     end
     'shoot' : begin
      child = widget_info(ev.top, find_by_uname='input_ra') 
      widget_control,child,get_value=stat,/no_copy
      status.ra=float(stat)

      child = widget_info(ev.top, find_by_uname='input_dec') 
      widget_control,child,get_value=stat,/no_copy
      status.dec=float(stat)

      ;; child = widget_info(ev.top, find_by_uname='input_oax') 
      ;; widget_control,child,get_value=stat,/no_copy
      ;; nu.oa[0]=float(stat)

      ;; child = widget_info(ev.top, find_by_uname='input_oay') 
      ;; widget_control,child,get_value=stat,/no_copy
      ;; nu.oa[1]=float(stat)

      if(sqrt(total((nu.oa-nu.oa_prev)^2)) gt 0.1) then begin
         ;; print,'Reset pointing for new OA:',nu.oa[0],nu.oa[1]
         ;; pos_prev=mm2det(nu.oa_prev)
         ;; pa = status.pa
         ;; make_astr_nustar, status.ra, status.dec, pa, astr=astr, oa=pos_prev
         ;; pos_new=mm2det(nu.oa)
         ;; xy2ad,pos_new[0],pos_new[1],astr,new_ra,new_dec
         ;; print,'New RA/Dec: ',new_ra,new_dec
         ;; status.ra=new_ra
         ;; status.dec=new_dec

         ;; child = widget_info(ev.top, find_by_uname='input_ra') 
         ;; widget_control,child,set_value=Number_Formatter(status.ra,Decimals=4),/no_copy
         
         ;; child = widget_info(ev.top, find_by_uname='input_dec') 
         ;; widget_control,child,set_value=Number_Formatter(status.dec,Decimals=4),/no_copy
      endif
      if(abs(nu.fov_shift_x) gt 0.0 or abs(nu.fov_shift_y) gt 0.0) then begin
         make_astr_nustar, status.ra, status.dec, status.pa, astr=astr_cur, oa=mm2det(nu.oa)
         xy2ad,mm2det(nu.oa[0]+nu.fov_shift_x)-1,mm2det(nu.oa[1]+nu.fov_shift_y)-1,astr_cur,new_ra,new_dec
         print,'New RA/Dec: ',new_ra,new_dec,' was:',status.ra, status.dec
         ;make_astr_nustar, new_ra, new_dec, status.pa, astr=astr, oa=mm2det(nu.oa)
         ;ad2xy,status.ra, status.dec, astr_cur, oax,oay
         ;print,'mm2det:',mm2det(nu.oa[0]),mm2det(nu.oa[1])
         ;print,'OA position and index',nu.oa[0],nu.oa[1],' ->',oax,oay

         status.ra=new_ra
         status.dec=new_dec
         child = widget_info(ev.top, find_by_uname='input_ra') 
         widget_control,child,set_value=cgNumber_Formatter(status.ra,Decimals=4),/no_copy
         
         child = widget_info(ev.top, find_by_uname='input_dec') 
         widget_control,child,set_value=cgNumber_Formatter(status.dec,Decimals=4),/no_copy
         nu.fov_shift_x = 0.0D 
         nu.fov_shift_y = 0.0D
      end

      stray_light_render

      nu.oa_prev=nu.oa
   end
   'up' : begin
      nu.fov_shift_y+=nu.fov_shift_step
      plot_targets
   end
   'upl' : begin
      nu.fov_shift_x-=nu.fov_shift_step
      nu.fov_shift_y+=nu.fov_shift_step
      plot_targets
   end
   'upr' : begin
      nu.fov_shift_x+=nu.fov_shift_step
      nu.fov_shift_y+=nu.fov_shift_step
      plot_targets
   end
   'downr' : begin
      nu.fov_shift_x+=nu.fov_shift_step
      nu.fov_shift_y-=nu.fov_shift_step
      plot_targets
   end
   'downl' : begin
      nu.fov_shift_x-=nu.fov_shift_step
      nu.fov_shift_y-=nu.fov_shift_step
      plot_targets
   end
   'down' : begin
      nu.fov_shift_y-=nu.fov_shift_step
      plot_targets
   end
   'left' : begin
      nu.fov_shift_x-=nu.fov_shift_step
      plot_targets
   end
   'right' : begin
      nu.fov_shift_x+=nu.fov_shift_step
      plot_targets
   end
   'center' : begin
      nu.fov_shift_x=0.0
      nu.fov_shift_y=0.0
      status.ra=status.ra_orig
      status.dec=status.dec_orig
      child = widget_info(ev.top, find_by_uname='input_ra') 
      widget_control,child,set_value=cgNumber_Formatter(status.ra,Decimals=4),/no_copy
      child = widget_info(ev.top, find_by_uname='input_dec') 
      widget_control,child,set_value=cgNumber_Formatter(status.dec,Decimals=4),/no_copy
      plot_targets
   end
   'null' : print,'Press "None"!'
   ;; 'oax' : print,'Press "Shoot"!'
   ;; 'oay' : print,'Press "Shoot"!'
   'ra' : print,'Press "Shoot"!'
   'dec' : print,'Press "Shoot"!'
   'draw' : ;widget_control, ev.top, /destroy
   'RequestPA' : print,'Request PA for ra/dec'
   'quit' : widget_control, ev.top, /destroy
 END 
END


;; smooth=0 -- no smoothing
;; smooth=1 -- take smoothing factors from stray_catalog
;; smooth>1 -- take this smooth factor


pro nuplan, pnt_ra, pnt_dec, pa=pnt_pa , stray_catalog=stray_catalog, oa=oa, oaa=oaa, oab=oab, smooth=smooth, $
            silent=silent, quit=quit, target_catalog=target_catalog, scan_step=scan_step, badpix=badpix, $
            key=key,slp_level=slp_level, fmin=fmin
common nuplan, nu, status, sources, target
forward_function read_combined_catalog

  if(n_elements(smooth) eq 0) then smooth=0  
  if(n_elements(key) eq 0) then key='nuplan'  
  if(n_elements(badpix) eq 0) then badpix=0  
  if(n_elements(oa) eq 0) then oa=[3., 3.] ; "default" position of OA in mm, will be converted to crpix WCS keyword
  if(n_elements(oaa) eq 0) then oaa=[3., 3.] 
  if(n_elements(oab) eq 0) then oab=[3., 3.] 
  if(n_elements(scan_step) eq 0) then scan_step=10.
  if(n_elements(silent) eq 0) then silent=0
  if(n_elements(quit) eq 0) then quit=0
  if n_elements(fmin) eq 0 then fmin = 5. ; default to 5 mCrab 


;npar = N_params()
;if ( npar eq 0 ) then begin
; print,'NuPLAN Usage: nuplan, ra, dec, pa, source_input_catalog'
; print,'ra,dec - pointing position, pa -- angle to the North Pole (PA).'          
; print, 'source_input_catalog is path to the source catalog in format: "name,ra,dec,flag"'
; print, 'All parameters can be changed later in program GUI.' 
; print, ''
;endif

  target={src_name:[''],src_ra:[0.0],src_dec:[0.0],src_flag:[0],index:[1]}

  if(n_elements(slp_level) eq 0) then slp_level=0      
  if(n_elements(pnt_ra) eq 0) then read,'NuPLAN: Enter R.A.: ',pnt_ra      
  if(n_elements(pnt_dec) eq 0) then read,'NuPLAN: Enter Dec.: ',pnt_dec      
  if(n_elements(pnt_pa) eq 0) then read,'NuPLAN: Enter PA: ',pnt_pa 
  ;if(n_elements(stray_catalog) eq 0) then read,'NuPLAN: Enter path to the source input catalog: ',stray_catalog 

  if(n_elements(target_catalog) gt 0) then begin
     read_source_catalog,target_catalog,target_name,target_ra,target_dec,target_flag,index=target_index
     target={src_name:target_name,src_ra:target_ra,src_dec:target_dec,src_flag:target_flag,index:target_index}
  endif 
  
hgap = 0.15 ; detector half gap [mm]                                                                                                          
n_detx = 64 ; number of detector x bins 
n_dety = 64 ; number of detector y bins 
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

print,'Saving [x,y]dim_mm.fits'
writefits, 'xdim_mm.fits',xpos_2d_array
writefits, 'ydim_mm.fits',ypos_2d_array


status={ra:pnt_ra, dec:pnt_dec, ra_orig:pnt_ra, dec_orig:pnt_dec, pa:pnt_pa,slpa:0,slpb:0,eff:0,vis:0,winid:0,mainid:0,infoid:0,TopTab:0,table:0,infotext:0,silent:silent,scan_step:scan_step, key:key, smooth:smooth}
nu={hgap:hgap,n_detx:n_detx,n_dety:n_dety,xpos_array:xpos_array,ypos_array:ypos_array,oa:oa,oaa:oaa,oab:oab,oa_prev:oa,fov_shift_x:0.0D,fov_shift_y:0.0D,fov_shift_step:0.5D,dr:!PI/180.,rd:180./!PI}

; Read in the combined BAT 70 Month and INTEGRAL galactic plane survey catalogs
read_combined_catalog, ra=pnt_ra, dec=pnt_dec, src_name, src_ra, src_dec, src_flux, src_flag, $
                fmin=fmin
 
;; if(n_elements(stray_catalog) gt 0) then begin
;;    read_source_catalog,stray_catalog,src_name,src_ra,src_dec,src_flag,index=index
;; endif
sources={src_name:src_name,src_ra:src_ra,src_dec:src_dec,src_flux:src_flux,src_flag:src_flag}

main = widget_base (title='NuPLAN', /column, xsize=1024, ysize=690)             ; main base
wTab = WIDGET_TAB(main, LOCATION=location)
status.TopTab=wTab
wT1 = WIDGET_BASE(wTab, TITLE='Stray light pattern', uvalue='tab1', /COLUMN) 

cntl = widget_base (wT1, /row, /frame)  
draw = widget_draw (wT1, uvalue='draw', /button, xsize=1024, ysize=500)               ; graphics pane
btn1 = widget_button (cntl, uvalue='quit', value='Quit')             ; quit button
wLabel_status = WIDGET_LABEL(cntl, uname='infotext', /ALIGN_LEFT, $
                             VALUE='                                                          ') 
status.infotext=wLabel_status

wLabel_info = WIDGET_LABEL(cntl, VALUE='') 
status.infoid=wLabel_info
;btn_request = widget_button (cntl, uvalue='RequestPA', value='Request PA')             ; quit button
cntl2=widget_base(wT1, /row, /frame)

cntl_pos=widget_base(cntl2, /column, /frame)
cntl_pos1=widget_base(cntl_pos, /row)
cntl_pos2=widget_base(cntl_pos, /row)

wLabel_pos = WIDGET_LABEL(cntl_pos1, VALUE="OA J2000 position (RA,Dec):") 
input_ra = widget_text (cntl_pos2, uname='input_ra', value=cgNumber_Formatter(status.ra,Decimals=4), uval='ra',/NO_NEWLINE,/EDITABLE, xsize=8) 
input_dec = widget_text (cntl_pos2, uname='input_dec', value=cgNumber_Formatter(status.dec,Decimals=4), uval='dec',/NO_NEWLINE,/EDITABLE, xsize=8) 

;; cntl_oa=widget_base(cntl2, /column, /frame)
;; cntl_oa1=widget_base(cntl_oa, /row)
;; cntl_oa2=widget_base(cntl_oa, /row)

;; wLabel_oa = WIDGET_LABEL(cntl_oa1, VALUE="OA det. position in mm:") 
;; input_ra = widget_text (cntl_oa2, uname='input_oax', value=Number_Formatter(nu.oa[0]), uval='oax',/NO_NEWLINE,/EDITABLE, xsize=6) 
;; input_dec = widget_text (cntl_oa2, uname='input_oay', value=Number_Formatter(nu.oa[1]), uval='oay',/NO_NEWLINE,/EDITABLE, xsize=6) 



cntl_fov=widget_base(cntl2, /column, /frame)
cntl_fov1=widget_base(cntl_fov, /row)
cntl_fov2=widget_base(cntl_fov, /row)
cntl_fov3=widget_base(cntl_fov, /row)

btn_up1 = widget_button (cntl_fov1, uvalue='upl',    value='        ')
btn_up  = widget_button (cntl_fov1, uvalue='up',     value='   Up   ')
btn_up2 = widget_button (cntl_fov1, uvalue='upr',    value='        ')

btn_md1 = widget_button (cntl_fov2, uvalue='left',   value='  Left  ')
btn_md  = widget_button (cntl_fov2, uvalue='center', value='  Reset ')
btn_md2 = widget_button (cntl_fov2, uvalue='right',  value='  Right ')

btn_bm1 = widget_button (cntl_fov3, uvalue='downl',  value='        ')
btn_bm  = widget_button (cntl_fov3, uvalue='down',   value='  Down  ')
btn_bm2 = widget_button (cntl_fov3, uvalue='downr',  value='        ')

cntl_pa=widget_base(cntl2, /column, /frame)
cntl_pa1=widget_base(cntl_pa, /row)
cntl_pa2=widget_base(cntl_pa, /row)

sld = widget_slider (cntl_pa2, min=0, max=360, value=status.pa, uval='pa', SCR_XSIZE=450)     ; slider 
btn2 = widget_button (cntl_pa2, uvalue='shoot', value='Shoot')             ; shoot button
btn3 = widget_button (cntl_pa2, uvalue='scan', value='Scan')             ; scan button
wLabel_sld = WIDGET_LABEL(cntl_pa1, VALUE="Positional Angle (PA)") 

wT2 = WIDGET_BASE(wTab, TITLE='Source catalog', uvalue='tab2', /COLUMN) 
  
  ;wSlider = WIDGET_SLIDER(wT2) 

  col_labels = ['Source', 'R.A.', 'Dec.', 'Flux (mCrab)','Flag'] 
  row_labels = [STRING(1,FORMAT='(i)')] 

  d={name:src_name[0],ra:src_ra[0],dec:src_dec[0],flux:src_flux[0],flag:src_flag[0]}
  data=[d]
  n_src=n_elements(src_name)
  for i=1,n_src-1 do begin
     d={source:src_name[i],ra:src_ra[i],dec:src_dec[i],flux:src_flux[i],flag:src_flag[i]}
     data=[data,d]
     row_labels=[row_labels,STRING(i+1,FORMAT='(i)')]
  endfor


  ; Combine structure data into a vector of structures. 
  ;data = [d1,d2]

  ; To make sure the table looks nice on all platforms, 
  ; set all column widths to the width of the longest string 
  ; that can be a header. 
  maxw=dblarr(4)
  max_strlen = strlen('  IGR J16283-4838  ') 
  maxwidth = max_strlen * !d.x_ch_size + 6   ; ... + 6 for padding 
  maxw[0]=maxwidth
  max_strlen = strlen('J16283-4838') 
  maxwidth = max_strlen * !d.x_ch_size + 6   ; ... + 6 for padding 
  maxw[1:3]=maxwidth
  
  table1 = WIDGET_TABLE(wT2, VALUE=data, /ROW_MAJOR, /ALL_EVENTS, $ 
    ROW_LABELS=row_labels, COLUMN_LABELS=col_labels, /editable, $ 
    COLUMN_WIDTHS=maxw, /RESIZEABLE_COLUMNS) 
  status.table=table1
wLabel = WIDGET_LABEL(wT2, VALUE="Don't forget to release cell after editing (e.g. by typing return button)") 
widget_control, main, /realize 
widget_control, draw, get_value=win_id ;get the window id
status.winid=win_id
status.mainid=main

stray_light_render, badpix=badpix

if(quit eq 1) then begin
   widget_control, main, /destroy
   return
endif

xmanager, 'nuplan', main, /no_block                         ; wait for events
end
