pro stray_light_render, badpix=badpix
  common nuplan, nu, status, sources, target

  if(n_elements(badpix) eq 0) then badpix=0    

  n_detx=nu.n_detx
  n_dety=nu.n_dety
  dr=nu.dr
  rd=nu.rd
  pa = status.pa
  hgap=nu.hgap

  RA_SRC_IN=sources.src_ra
  DEC_SRC_IN=sources.src_dec
  
  RA_SRC = RA_SRC_IN*dr
  DEC_SRC = DEC_SRC_IN*dr
 

 ; (BG)  Assumes that src_flux is now in some reasonable units
;  Flux = replicate(1., n_elements(RA_SRC))
  flux = sources.src_flux

  n_src = n_elements(RA_SRC) 

   ; Kaya's code. Make mask
   dmask = fltarr(n_detx+n_dety, n_detx)
   dmask0 = fltarr(n_detx+n_dety, n_detx)
   ; cycle on sources

   dmask_fp1=fltarr(n_detx, n_dety) 
   dmask_fp2=fltarr(n_detx, n_dety) 
   
   file_delete,'slp_'+status.key+'_fpA_raw.reg',/ALLOW_NONEXISTENT
   file_delete,'slp_'+status.key+'_fpB_raw.reg',/ALLOW_NONEXISTENT
   file_delete,'slp_'+status.key+'_fpA_sky.reg',/ALLOW_NONEXISTENT
   file_delete,'slp_'+status.key+'_fpB_sky.reg',/ALLOW_NONEXISTENT

   if (status.silent ne 1) then  print, 'i, name, RA [deg], DEC [deg], OAA [deg], AZ [deg], Flux (mCrab) =
   for ii=0, n_src-1 do begin 
      if(sources.src_flag[ii] eq 0) then continue
      OAA_var = arclength(status.ra*dr, status.dec*dr, RA_SRC(ii), DEC_SRC(ii))
      az_angle = AZIMUTH_ANGLE(RA_SRC(ii), DEC_SRC(ii), status.ra*dr, status.dec*dr)
; (BG) Added some more details that are useful for planning:
;      if(status.silent ne 1) then print, 'i, RA [deg], DEC [deg], OAA [deg], AZ [deg] = ', ii+1, RA_SRC(ii)*rd, DEC_SRC(ii)*rd, OAA_var*rd, az_angle*rd
      if(status.silent ne 1) then print,string(ii+1, format ='(i0)'),' ', strtrim(sources.src_name[ii]), $
                                         string(RA_SRC(ii)*rd, format = '(d10.3)'), DEC_SRC(ii)*rd, OAA_var*rd, az_angle*rd, ' ', sources.src_flux[ii]

      dmask0=LEAKAGE_MAP(OAA_var, az_angle, pa*dr)*Flux(ii)

      ;if(status.silent ne 1) then begin
      if(1) then begin
   ; find index of the OA given in mm
         ;get_oa,nu.oa, pos=oa_index
         make_astr_nustar, status.ra, status.dec, pa, astr=astr, oa=mm2det(nu.oa)

         make_astr_nustar, status.ra, status.dec, pa, astr=astra, oa=mm2det(nu.oaa)
         make_astr_nustar, status.ra, status.dec, pa, astr=astrb, oa=mm2det(nu.oab)

         dmask0_fp1=fltarr(n_detx, n_dety) 
         dmask0_fp2=fltarr(n_detx, n_dety) 
                                ; Kaya's code.
         for i=0, n_detx-1 do begin 
            for j=0, n_dety-1 do begin
               dmask0_fp1(i,j) = (dmask0(i,j) gt 0.0)? dmask0(i, j) : 0.0
               dmask0_fp2(i,j) = (dmask0(i+n_detx,j) gt 0.0)? dmask0(i+n_detx, j) : 0.0
            endfor 
         endfor

         if (sources.src_flag[ii] ge 1 and status.smooth ge 1) then begin
           
            smooth_factor=sources.src_flag[ii]
            if(status.smooth gt 1) then smooth_factor=status.smooth

            dmask0_fp1s=smooth(dmask0_fp1,smooth_factor,/EDGE_TRUNCATE)
            index=where(dmask0_fp1s gt 0.01, cnt)
            if(cnt gt 0) then dmask0_fp1(index)=1

            dmask0_fp2s=smooth(dmask0_fp2,smooth_factor,/EDGE_TRUNCATE)
            index=where(dmask0_fp2s gt 0.01, cnt)
            if(cnt gt 0) then dmask0_fp2(index)=1

            writefits, 'dmask_fpA_sm.fits',dmask0_fp1s,hdr
            writefits, 'dmask_fpB_sm.fits',dmask0_fp2s,hdr
         endif else begin
         endelse
 
         Contour, dmask0_fp1, PATH_INFO=info, PATH_XY=xy, XSTYLE=1, YSTYLE=1, /PATH_DATA_COORDS, /CLOSED, /NLEVELS
         if(n_elements(xy) gt 0) then begin
         openw,1,'slp_'+status.key+'_fpA_raw.reg',/APPEND
         printf,1,'image; polygon(',format='(a, $)' 
         openw,2,'slp_'+status.key+'_fpA_sky.reg',/APPEND
         printf,2,'fk5; polygon(',format='(a, $)' 
         for t=0,n_elements(xy[0,*])-1 do begin
            xpix=xy[0,t]+1.
            ypix=xy[1,t]+1.
            printf,1,xpix,',',ypix,format='(f12.4,a,f12.4, $)' 
            xy2ad,xpix-1.,ypix-1.,astra,ra,dec
            printf,2,ra,',',dec,format='(f12.4,a,f12.4,a, $)' 
         endfor
         printf,1,') # color=cyan width=3 text={'+sources.src_name(ii)+'}'
         printf,2,') # color=cyan width=3 text={'+sources.src_name(ii)+'}'
         close,1
         close,2
; (BG) Construct labels for plotting below
         label_x = mean(xy[0, *])
         label_y = mean(xy[1, *])
         this_src_name = sources.src_name[ii]
         push, labels1_x, label_x
         push, labels1_y, label_y
         push, labels1_name, this_src_name


         endif

         Contour, dmask0_fp2, PATH_INFO=info, PATH_XY=xy, XSTYLE=1, YSTYLE=1, /PATH_DATA_COORDS, /CLOSED, /NLEVELS
         if(n_elements(xy) gt 0) then begin
         openw,1,'slp_'+status.key+'_fpB_raw.reg',/APPEND
         printf,1,'image; polygon(',format='(a, $)' 
         openw,2,'slp_'+status.key+'_fpB_sky.reg',/APPEND
         printf,2,'fk5; polygon(',format='(a, $)' 
         for t=0,n_elements(xy[0,*])-1 do begin
            xpix=xy[0,t]+1.
            ypix=xy[1,t]+1.
            printf,1,xpix,',',ypix,format='(f12.4,a,f12.4, $)' 
            xy2ad,xpix-1.,ypix-1.,astrb,ra,dec
            printf,2,ra,',',dec,format='(f12.4,a,f12.4,a, $)' 
         endfor
         printf,1,') # color=cyan width=3 text={'+sources.src_name(ii)+'}'
         printf,2,') # color=cyan width=3 text={'+sources.src_name(ii)+'}'
         close,1
         close,2

; (BG) Dittor for FPMB:
         label_x = mean(xy[0, *])
         label_y = mean(xy[1, *])
         this_src_name = sources.src_name[ii]
         push, labels2_x, label_x
         push, labels2_y, label_y
         push, labels2_name, this_src_name


         endif
      endif

   dmask_fp1+=dmask0_fp1 
   dmask_fp2+=dmask0_fp2 
   ;dmask+=dmask0 old
   endfor
   
; old:
;   for i=0, n_detx-1 do begin 
;      for j=0, n_dety-1 do begin 
;         dmask_fp1(i,j) = (dmask(i,j) gt 0.0)? dmask(i, j) : 0.0
;         dmask_fp2(i,j) = (dmask(i+n_detx,j) gt 0.0)? dmask(i+n_detx, j) : 0.0
;      endfor 
;   endfor 

if(1) then begin 
   ; find index of the OA given in mm
   ;get_oa,nu.oa, pos=oa_index
   make_header_nustar, status.ra, status.dec, pa, astr=astr, hdr=hdr, bitpix=4, oa=mm2det(nu.oa)
   make_header_nustar, status.ra, status.dec, pa, astr=astra, hdr=hdra, bitpix=4, oa=mm2det(nu.oaa)
   make_header_nustar, status.ra, status.dec, pa, astr=astrb, hdr=hdrb, bitpix=4, oa=mm2det(nu.oab)
   ;ad2xy,status.ra, status.dec, astr, oax,oay
   ;print,'mm2det:',mm2det(nu.oa[0]),mm2det(nu.oa[1])
   ;print,'OA position and index',nu.oa[0],nu.oa[1],' ->',oax,oay
   print,'Saving SLP sources to slp_sources.reg'
   openw, lun,  'slp_sources.reg', /get_lun
   for ii=0, n_src-1 do begin 
      if(sources.src_flag[ii] eq 0) then continue
      printf,lun,'fk5; circle( ',$
             sources.src_ra(ii),', ',sources.src_dec(ii),', 0.15 ) # text={'+sources.src_name(ii)+'}'
   endfor
   close, lun
   free_lun, lun

   ; check is there hotpix matrix
   filesav='stack_hotpix.sav'
   if(file_test(filesav)) then begin
      RESTORE, filesav
      dmask_fp1+=hotpix_raw[*,*,0]
      dmask_fp2+=hotpix_raw[*,*,1]
   endif

   ; check is there badpix-region matrix
   filesav='region_badpix.sav'
   if(file_test(filesav)) then begin
      print,'stray_light_render: RESTORE ',filesav
      RESTORE, filesav
      dmask_fp1+=hotpix_raw[*,*,0]
      dmask_fp2+=hotpix_raw[*,*,1]
      writefits,'masked_RAW_tot.fpA.img',dmask_fp1
      writefits,'masked_RAW_tot.fpB.img',dmask_fp2
   endif

   print,'Saving dmask_fp[A,B].fits'
   dmask_fp1[5,34]+=1
   dmask_fp1[5,35]+=1
   dmask_fp1[5,33]+=1
   dmask_fp1[4,34]+=1
   dmask_fp1[6,34]+=1
   writefits, 'dmask_fpA.fits',dmask_fp1,hdra
   writefits, 'dmask_fpB.fits',dmask_fp2,hdrb
      
   ; make user-defined badpixel matrix
   if(badpix eq 1) then begin
      print,"Run badpix code. Is Brian's code ready to use?"
      make_badpix,dmask_fp1,dmask_fp2
      ;; Stupid, but what we can do? The static name of this file
      ;; emerges from Brian's script.
      FILE_MOVE, 'nuAuserbadpix20100101v002.fits', 'nuAuserbadpix20100101v002.'+status.key+'.fits', /NOEXPAND_PATH, /OVERWRITE
      FILE_MOVE, 'nuBuserbadpix20100101v002.fits', 'nuBuserbadpix20100101v002.'+status.key+'.fits', /NOEXPAND_PATH, /OVERWRITE
   endif
endif

index=where(dmask(0:63,*) gt 0.0,count)
fp1_pct=0.0
if(count gt 0) then fp1_pct=count*100.0/64.^2

index=where(dmask(64:127,*) gt 0.0,count)
fp2_pct=0.0
if(count gt 0) then fp2_pct=count*100.0/64.^2
 
   index=where(dmask(0:31,0:31) gt 0.0,count)
   fp1chip2_pct=0.0
   if(count gt 0) then fp1chip2_pct=count*100.0/32.^2
   index=where(dmask(0:31,32:63) gt 0.0,count)
   fp1chip1_pct=0.0
   if(count gt 0) then fp1chip1_pct=count*100.0/32.^2
   index=where(dmask(32:63,32:63) gt 0.0,count)
   fp1chip0_pct=0.0
   if(count gt 0) then fp1chip0_pct=count*100.0/32.^2
   index=where(dmask(32:63,0:31) gt 0.0,count)
   fp1chip3_pct=0.0
   if(count gt 0) then fp1chip3_pct=count*100.0/32.^2

   index=where(dmask(64:95,0:31) gt 0.0,count)
   fp2chip2_pct=0.0
   if(count gt 0) then fp2chip2_pct=count*100.0/32.^2
   index=where(dmask(64:95,32:63) gt 0.0,count)
   fp2chip1_pct=0.0
   if(count gt 0) then fp2chip1_pct=count*100.0/32.^2
   index=where(dmask(96:127,32:63) gt 0.0,count)
   fp2chip0_pct=0.0
   if(count gt 0) then fp2chip0_pct=count*100.0/32.^2
   index=where(dmask(96:127,0:31) gt 0.0,count)
   fp2chip3_pct=0.0
   if(count gt 0) then fp2chip3_pct=count*100.0/32.^2

   fp1or2_pct=0.0
   index=where (dmask_fp1 eq 0.0 or dmask_fp2 eq 0.0, count)
   if(count gt 0) then fp1or2_pct=count*100.0/64.^2
   if(status.silent ne 1) then $
      print, 'FPA, FPB, 200-FPA-FPB, FPA.or.FPB [%]= ', fp1_pct, fp2_pct, (100-fp1_pct)+(100-fp2_pct), fp1or2_pct,format='(a,4f8.2)'
   vis=0.0
   index=where((dmask_fp1*dmask_fp2) eq 0, count)
   if(count gt 0) then vis=count*100.0/64.^2
   eff=(100-fp1_pct)+(100-fp2_pct)
   Label = widget_info(status.mainid, find_by_uname='infotext') 
   widget_control,Label,set_value= $
                  'Efficiency: '+cgNumber_Formatter(eff,Decimals=1)+ $
                  '% Visibility: '+cgNumber_Formatter(vis,Decimals=1)+ '%'
   status.vis=vis
   status.eff=eff
   status.slpa=fp1_pct
   status.slpb=fp2_pct

   !p.multi=[0,2,1]
   !p.charsize=1.8
   bb=20.0
   contour, dmask_fp1, nu.xpos_array, nu.ypos_array, /cell_fill, xtit='DETX [mm]', ytit='DETY [mm]', $
            tit='FPA ', $
            levels=sources.src_flux[sort(sources.src_flux)] ; (BG) Contour levesl == source fluxes.
            ;,levels=[0,1,2,3,4]
            ;tit='FPA '+'(SLP '+string(fp1_pct, format='(f4.0)')+'%)',levels=[0,1,2,3,4]

   oplot, [nu.oa[0]], [nu.oa[1]], psym=2    
   oplot, [hgap,hgap], [-bb,bb], linestyle=0
   oplot, [-hgap,-hgap], [-bb,bb], linestyle=0
   oplot, [-bb,bb], [hgap,hgap], linestyle=0
   oplot, [-bb,bb], [-hgap,-hgap], linestyle=0
   ;xyouts, 1. , 17., 'chip0: '+string(fp1chip0_pct, format='(f4.0)')+'%'
   ;xyouts, 1. , -3., 'chip3: '+string(fp1chip3_pct, format='(f4.0)')+'%'
   ;xyouts, -18. , -3., 'chip2: '+string(fp1chip2_pct, format='(f4.0)')+'%'
   ;xyouts, -18. , 17., 'chip1: '+string(fp1chip1_pct, format='(f4.0)')+'%'
; (BG) Add labels:
   if n_elements(labels1_name) gt 0 then begin
      for ll = 0, n_elements(labels1_x)-1 do begin
         ; Convert to mm:
         mm_x = nu.xpos_array[round(labels1_x[ll])]
         mm_y = nu.ypos_array[round(labels1_y[ll])]
         xyouts, mm_x, mm_y, labels1_name[ll], color = cgColor('Cyan'), /data
      endfor
   endif



   green = GETCOLOR('green', 100)
   for jj=0, n_elements(target.src_ra)-1 do begin 
      if(target.src_flag[jj] eq 0) then continue
      ad2xy,target.src_ra[jj],target.src_dec[jj],astr,dx,dy
      if(dx gt 0 and dx lt n_detx and dy gt 0 and dy lt n_dety) then begin
         ;print,'TARGET:',target.src_name[jj],dx,dy,target.src_ra[jj],target.src_dec[jj],det2mm(dx),det2mm(dy)
         offset=sphdist(target.src_ra[jj],target.src_dec[jj],status.ra,status.dec,/DEGREES)*60.
         print,'TARGET OA offset, arcmin: ',target.src_name[jj],offset,format='(a,1x,a,1x,f8.2)'
      endif else begin
         continue
      endelse
      push,ptx,det2mm(dx+1)
      push,pty,det2mm(dy+1)
      push,ptn,target.src_name[jj]
      push,off,offset
   endfor
   
   for n=0,n_elements(ptx)-1 do begin
      PLOTS, CIRCLE(ptx[n], pty[n], 2.5), color=green
      oplot, [ptx[n]], [pty[n]], psym=1, color=green
      xyouts, ptx[n] , pty[n]+3.0, ptn[n], CHARSIZE=0.9, color=green
      xyouts, ptx[n]+3.0 , pty[n], String(off[n],format='(f3.1)')+"'", CHARSIZE=0.9, color=green
   endfor

   contour, dmask_fp2, nu.xpos_array, nu.ypos_array, /cell_fill, xtit='DETX [mm]', ytit='DETY [mm]', $
            tit='FPB ', $
            levels=sources.src_flux[sort(sources.src_flux)] ; (BG) Contour levesl == source fluxes.
            ;levels=[0,1,2,3,4]
            ;tit='FPB '+'(SLP '+string(fp2_pct,format='(f4.0)')+'%)',levels=[0,1,2,3,4]

   oplot, [nu.oa[0]], [nu.oa[1]], psym=2
   oplot, [hgap,hgap], [-bb,bb], linestyle=0
   oplot, [-hgap,-hgap], [-bb,bb], linestyle=0
   oplot, [-bb,bb], [hgap,hgap], linestyle=0
   oplot, [-bb,bb], [-hgap,-hgap], linestyle=0
   ;xyouts, 1. , 17., 'chip0: '+string(fp2chip0_pct, format='(f4.0)')+'%'
   ;xyouts, 1. , -3., 'chip3: '+string(fp2chip3_pct, format='(f4.0)')+'%'
   ;xyouts, -18. , -3., 'chip2: '+string(fp2chip2_pct, format='(f4.0)')+'%'
   ;xyouts, -18. , 17., 'chip1: '+string(fp2chip1_pct, format='(f4.0)')+'%'

; (BG) Add labels:
   if n_elements(labels2_name) gt 0 then begin
      for ll = 0, n_elements(labels2_x)-1 do begin
         ; Convert to mm:
         mm_x = nu.xpos_array[round(labels2_x[ll])]
         mm_y = nu.ypos_array[round(labels2_y[ll])]
         xyouts, mm_x, mm_y, labels2_name[ll], color = cgColor('Cyan'), /data
      endfor
   endif


   for n=0,n_elements(ptx)-1 do begin
      PLOTS, CIRCLE(ptx[n], pty[n], 2.5), color=green
      oplot, [ptx[n]], [pty[n]], psym=1, color=green
      xyouts, ptx[n] , pty[n]+3.0, ptn[n], CHARSIZE=0.9, color=green
      xyouts, ptx[n]+3.0 , pty[n], String(off[n],format='(f3.1)')+"'", CHARSIZE=0.9, color=green
   endfor

end
