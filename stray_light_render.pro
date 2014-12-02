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
  
                                ; (BG) Removed as out of date.
 
  if (status.silent ne 1) then  print, 'i, name, RA [deg], DEC [deg], OAA [deg], AZ [deg], Flux (mCrab) ='
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
        ;; make_astr_nustar, status.ra, status.dec, pa, astr=astr, oa=mm2det(nu.oa)

        ;; make_astr_nustar, status.ra, status.dec, pa, astr=astra, oa=mm2det(nu.oaa)
        ;; make_astr_nustar, status.ra, status.dec, pa, astr=astrb, oa=mm2det(nu.oab)

        dmask0_fp1=fltarr(n_detx, n_dety) 
        dmask0_fp2=fltarr(n_detx, n_dety) 
                                ; Kaya's code.
        ;; for i=0, n_detx-1 do begin 
        ;;    for j=0, n_dety-1 do begin
        ;;       dmask0_fp1(i,j) = (dmask0(i,j) gt 0.0)? dmask0(i, j) : 0.0
        ;;       dmask0_fp2(i,j) = (dmask0(i+n_detx,j) gt 0.0)? dmask0(i+n_detx, j) : 0.0
        ;;    endfor 
        ;; endfor
        dmask0_fp1 = dmask0[0:63, *]
        dmask0_fp2 = dmask0[64:*, *]



        
        Contour, dmask0_fp1, PATH_INFO=info, PATH_XY=xy, XSTYLE=1, YSTYLE=1, /PATH_DATA_COORDS, /CLOSED, /NLEVELS
        if(n_elements(xy) gt 0) then begin
 
; (BG) Construct labels for plotting below
           label_x = mean(xy[0, *])
           label_y = mean(xy[1, *])
           this_src_name = sources.src_name[ii]
           push, flux1, sources.src_flux[ii]
           push, labels1_x, label_x
           push, labels1_y, label_y
           push, labels1_name, this_src_name
           
        endif

        undefine, xy
        Contour, dmask0_fp2, PATH_INFO=info, PATH_XY=xy, XSTYLE=1, YSTYLE=1, /PATH_DATA_COORDS, /CLOSED, /NLEVELS
        if(total(dmask0_fp2) gt 0) then begin

; (BG) Ditto for FPMB:
           label_x = mean(xy[0, *])
           label_y = mean(xy[1, *])
           this_src_name = sources.src_name[ii]
           push, flux2, sources.src_flux[ii]
           push, labels2_x, label_x
           push, labels2_y, label_y
           push, labels2_name, this_src_name

           
           ;; index=where(dmask0_fp2 gt 0.0,count)
           ;; print, '' 
           ;; print, this_src_name, count * 100.0 / 64.^2
           ;; print, ''
           ;; stop


        endif
     endif

     dmask_fp1+=dmask0_fp1 
     dmask_fp2+=dmask0_fp2 
                                ;dmask+=dmask0 old
  endfor



  index=where(dmask_fp1 gt 0.0,count)
  fp1_pct=0.0
  if(count gt 0) then fp1_pct=count*100.0/64.^2

  index=where(dmask_fp2 gt 0.0,count)
  fp2_pct=0.0
  if(count gt 0) then fp2_pct=count*100.0/64.^2
  

  index=where(dmask_fp1[32:*, 32:*] gt 0.0,count)
  fp1chip0_pct=0.0
  if(count gt 0) then fp1chip0_pct=count*100.0/32.^2
  index=where(dmask_fp2[32:*, 32:*] gt 0.0,count)
  fp2chip0_pct=0.0
  if(count gt 0) then fp2chip0_pct=count*100.0/32.^2

  if(status.silent ne 1) then begin
     print, 'Det0 Loss (%), FPMA, FPMB: ', fp1chip0_pct, fp2chip0_pct
     print, 'All Dets Loss (%): FPA, FPB:', fp1_pct, fp2_pct
     print, 'Average available (%):', 0.5*((100-fp1_pct)+(100-fp2_pct)), format='(a,4f8.2)'
  endif 

  

  vis=0.0
  index=where((dmask_fp1*dmask_fp2) eq 0, count)
  if(count gt 0) then vis=count*100.0/64.^2

  eff=(100-fp1_pct)+(100-fp2_pct)

  status.eff=eff
  status.slpa=fp1_pct
  status.slpb=fp2_pct
  status.loss0=fp1chip0_pct
  status.loss1=fp2chip0_pct

  !p.multi=[0,2,1]
  !p.charsize=1.25
  !p.color = cgcolor('black')
  !p.background = cgColor('White')
  bb=20.0
  sources.src_flux = sources.src_flux+randomu(seed, n_elements(sources.src_flux))* 0.001
  position = [0.1, 0.1, 0.4, 0.5]
;  cgloadct, 0, ncolors = n_elements(sources.src_flux)+1, bottom = 0, /reverse
;  c_colors = indgen(n_elements(sources.src_flux)+1)

  if max(dmask_fp1) gt 0 then begin

     
     cgloadct, 0, /reverse
     cgloadct, 0, ncolors = n_elements(flux1+1), clip=[50, 200], /reverse, bottom = 1
     c_colors = indgen(n_elements(flux1)+1)
     flux1 += randomu(seed, n_elements(flux1))*0.001


     flux_levels = flux1[sort(flux1)]

     cgcontour, dmask_fp1, nu.xpos_array, nu.ypos_array, /cell_fill, xtit='DETX [mm]', ytit='DETY [mm]', $
                tit='FPMA ', position=position,  /closed, c_colors = c_colors, levels = [0, flux_levels-0.01]


  endif else plot, nu.xpos_array, nu.ypos_array, /nodata, xtit='DETX [mm]', ytit='DETY [mm]', $
                    tit='FPMA ', position=position
  



;$
;           levels=[-1, sources.src_flux[sort(sources.src_flux)], max(sources.src_flux)*1.5] ; (BG) Contour levesl == source fluxes.
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
                                ;xyouts, -18. , 17., 'chip1:
                                ;'+string(fp1chip1_pct,
                                ;format='(f4.0)')+'%'


; (BG) Add labels:
  if n_elements(labels1_name) gt 0 then begin
     for ll = 0, n_elements(labels1_x)-1 do begin
                                ; Convert to mm:
        mm_x = nu.xpos_array[round(labels1_x[ll])]
        mm_y = nu.ypos_array[round(labels1_y[ll])]
       
        xyouts, mm_x, mm_y, labels1_name[ll], color = cgColor('Red'), /data
     endfor
  endif
  xyouts, /data, 21, 21, 'PA Angle '+string(pa, format = '(i0)')




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


  position = [0.6, 0.1, 0.9, 0.5]

  if max(dmask_fp2) gt 0 then begin
     cgloadct, 0, /reverse
     cgloadct, 0, ncolors = n_elements(flux2), clip=[50, 200], /reverse, bottom = 1
     c_colors = indgen(n_elements(flux2)+1)
     flux2 += randomu(seed, n_elements(flux2))*0.001

     flux_levels = flux2[sort(flux2)]

     cgcontour, dmask_fp2, nu.xpos_array, nu.ypos_array, /cell_fill, xtit='DETX [mm]', ytit='DETY [mm]', $
                tit='FPMB ', position=position,  /closed, c_colors = c_colors, levels = [0, flux_levels-0.01]



  endif else plot, nu.xpos_array, nu.ypos_array, /nodata, xtit='DETX [mm]', ytit='DETY [mm]', $
                    tit='FPMB ', position=position




;  contour, /iso, dmask_fp2, nu.xpos_array, nu.ypos_array, /cell_fill, xtit='DETX [mm]', ytit='DETY [mm]', $
;           tit='FPB ', $
;           levels=sources.src_flux[sort(sources.src_flux)] ; (BG) Contour levesl == source fluxes.
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
        xyouts, mm_x, mm_y, labels2_name[ll], color = cgColor('Red'), /data
     endfor
  endif


  for n=0,n_elements(ptx)-1 do begin
     PLOTS, CIRCLE(ptx[n], pty[n], 2.5), color=green
     oplot, [ptx[n]], [pty[n]], psym=1, color=green
     xyouts, ptx[n] , pty[n]+3.0, ptn[n], CHARSIZE=0.9, color=green
     xyouts, ptx[n]+3.0 , pty[n], String(off[n],format='(f3.1)')+"'", CHARSIZE=0.9, color=green
  endfor


end
