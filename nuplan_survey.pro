pro get_oa,oa,pos=pos
  common nuplan, nu

   if(nu.oa[0] gt max(nu.xpos_array) or $
      nu.oa[0] lt min(nu.xpos_array) or $
      nu.oa[1] gt max(nu.ypos_array) or $
      nu.oa[1] lt min(nu.ypos_array)) then begin
      print,'Error: OA is out of FOV'
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


pro nuplan_survey, pnt_file=pnt_file, map_center=map_center, stray_catalog=stray_catalog, silent=silent, key=key, pa_array=pa_array, oa=oa, slp_level=slp_level
Compile_Opt defint32

common nuplan, nu

  if(n_elements(slp_level) eq 0) then slp_level=3      

  if(n_elements(oa) eq 0) then oa=[3., 3.] ; "default" position of OA in mm, will be converted to crpix WCS keyword

; this is for personal use. If you don't have DATA_ROOT, don't worry
root = GETENV('DATA_ROOT')
; output directory (should exist)
saveto='SLP/'
if file_test(saveto) eq 0 then begin
FILE_MKDIR,saveto
endif

if(n_elements(silent) eq 0) then silent=0

pi=3.14159D
sq2pi=sqrt(2.0D*pi)
sq2=sqrt(2.0D)

; GRID1. Simply grab WCS from FITS file
;chandra_expomap = READFITS( root+'/Data/Norma/Chandra/expmap_hiE_2arcsec.fits', h, /NOSCALE)
;print,'use key: ',key
;sz=SIZE(chandra_expomap)
;wx=sz(1)
;wy=sz(2)

; read catalog
Deg2Rad = !PI/180.
Rad2Deg = 180./!PI 

; GRID2
; initiate detector grid (from Kaya's code)
; bitmask          
hgap = 0.15 ; detector half gap [mm]                                                                                                          
n_detx = 64 ; number of detector x bins 
n_dety = 64 ; number of detector y bins 
d0mask = fltarr(n_detx, n_dety)
d1mask = fltarr(n_detx, n_dety)
xpos_array = fltarr(n_detx)
ypos_array = fltarr(n_dety)
for xx = 0, n_detx-1 do begin
   for yy = 0, n_dety-1 do begin
      q_x = sign(xx+1 - 32.5)                                        
      q_y = sign(yy+1 - 32.5)
      xpos_array(xx) = (xx+1-32.5)*0.6048 + hgap*sign(xx+1-32.5)
      ypos_array(yy) = (yy+1-32.5)*0.6048 + hgap*sign(yy+1-32.5)
   endfor
endfor

nu={hgap:hgap,n_detx:n_detx,n_dety:n_dety,xpos_array:xpos_array,ypos_array:ypos_array,oa:oa,oa_prev:oa,dr:!PI/180.,rd:180./!PI}


; GRID3
; this is grid in galactic latitude
; this grid is used only for high-level products like SLP coverage in %
gal=map_center
gwx=20.5d; 2 degrees wide in both directions
gwy=20.5d
hnwx=1000; how many pixels in half of width
hnwy=1000
nwx=hnwx*2+1; the full number of pixels (this was done to have center of the map in the central pixel) check maps carefully in ds9!
nwy=hnwy*2+1
stepx=gwx/(nwx-1)
stepy=gwy/(nwy-1)
nn=nwx*nwy

map_fp1=dblarr(nwx,nwy); the map for FP1
map_fp2=dblarr(nwx,nwy); the map for FP2
map_fp1_exp=dblarr(nwx,nwy); the map for FP1
map_fp2_exp=dblarr(nwx,nwy); the map for FP2
map_cnt_fov_fpAB=lonarr(nwx,nwy); exposure map
map_cnt_fov_fpA=lonarr(nwx,nwy); exposure map
map_cnt_fov_fpB=lonarr(nwx,nwy); exposure map
map_cnt_pix=lonarr(nwx,nwy); exposure map

; make WCS coordinates for GRID3 (Note Cartesian projection -CAR)
  crpix=[fix(nwx/2.0d)+1,fix(nwy/2.0d)+1]
  cdelt=[-stepx,stepy]
  mkhdr, grid_hdr, 4, [nwx,nwy]
  make_astr,astr, DELTA = cdelt, CRPIX = crpix, $
            CRVAL = [gal(0),gal(1)], CTYPE=['GLON-CAR','GLAT-CAR'];, RADECSYS = 'FK5', EQUINOX = 2000.0
  putast, grid_hdr, astr, CD_TYPE=2

  read_nustar_pointings,pnt_file,pnt_ra,pnt_dec,pnt_lon=pnt_lon,pnt_lat=pnt_lat,obsid=obsid
  n_pnt=n_elements(pnt_ra)
  
  koeff=0.0d                    ; not used
  !p.multi=[0,2,1]

if (slp_level ge 0) then begin
   EULER, map_center[0], map_center[1], ra, dec, 2
   read_integral9, ra=ra, dec=dec, src_name, src_ra, src_dec, src_flag, $
                   slp_level=slp_level,  Rmax=22.0, Rmin=0.0
   index=where(src_flag eq 1)
endif


; read catalog of sources
if(n_elements(stray_catalog) gt 0) then begin
   read_source_catalog,stray_catalog,src_name,src_ra,src_dec,src_flag,index=index
endif

; save source list used in this run
openw,1,saveto+'sources.'+key+'.dat'
for ii=0, n_elements(src_name)-1 do if (src_flag(ii) eq 1) then printf,1,src_name(ii),',',src_ra(ii),', ',src_dec(ii),',',src_flag(ii)
close,1


;for l=0,n_elements(src_name)-1 do begin
;key=STRCOMPRESS(src_name[l],/REMOVE_ALL)
;RA_SRC_IN=[src_ra[l]]
;DEC_SRC_IN=[src_dec[l]]
;print,n_elements(index)
;if(n_src gt 0) then begin

RA_SRC_IN=src_ra(index)
DEC_SRC_IN=src_dec(index)

RA_SRC = RA_SRC_IN*Deg2Rad
DEC_SRC = DEC_SRC_IN*Deg2Rad
Flux = replicate(1., n_elements(RA_SRC))
n_src = n_elements(RA_SRC) 
;endif

file_delete,saveto+'PA.'+key+'.dat',/ALLOW_NONEXISTENT
; set PA array:

for r=0,n_elements(pa_array)-1 do begin
   PAdeg=pa_array[r]
   PA=PAdeg*Deg2Rad             ; in radians
   print,'PA ',PAdeg
   

   map_cnt_fov_fpAB(*,*)=0
   map_cnt_fov_fpA(*,*)=0
   map_cnt_fov_fpB(*,*)=0
   map_cnt_pix(*,*)=0
   map_fp1(*,*)=0.
   map_fp2(*,*)=0.

; running over survey individual pointings
for t=0,n_pnt-1 do begin
   ; take fixed PA angle
   ;PAdeg=2.0

      this_pnt_ra=pnt_ra(t)
      this_pnt_dec=pnt_dec(t)
   

   ; make header for individual observation
      get_oa, oa, pos=oa_index
      make_astr_nustar, this_pnt_ra, this_pnt_dec, PAdeg, astr=astr, oa=oa_index
      mkhdr, hdr_long, 3, [64,64]  
      putast, hdr_long, astr, CD_TYPE=2
      mkhdr, hdr_float, 4, [64,64]  
      putast, hdr_float, astr, CD_TYPE=2

      ;xy2ad,oa_index[0]-1,oa_index[1]-1,astr,ra, dec
      get_oa, [3.0,3.0], pos=oa_index_33
      xy2ad,oa_index_33[0]-1,oa_index_33[1]-1,astr,ra, dec
      print,obsid[t], this_pnt_ra, this_pnt_dec, ra, dec, format='(a,4f18.7)'

      ;continue

   ; Kaya's code. Make mask
   dmask = fltarr(n_detx+n_dety, n_detx)
   ; cycle on sources
   for ii=0, n_src-1 do begin 
      OAA = arclength(this_pnt_ra*Deg2Rad, this_pnt_dec*Deg2Rad, RA_SRC(ii), DEC_SRC(ii))
      az_angle = AZIMUTH_ANGLE(RA_SRC(ii), DEC_SRC(ii), this_pnt_ra*Deg2Rad, this_pnt_dec*Deg2Rad)
      if(silent ne 1) then print, 'i, RA [deg], DEC [deg], OAA [deg], AZ [deg] = ', $
                                  ii, RA_SRC(ii)*Rad2Deg, DEC_SRC(ii)*Rad2Deg, OAA*Rad2Deg, az_angle*Rad2Deg
      dmask+=LEAKAGE_MAP(OAA, az_angle, PA)*Flux(ii)
   endfor
   dmask_fp1=fltarr(n_detx, n_dety) 
   dmask_fp2=fltarr(n_detx, n_dety) 
   dmask_fp1_exp=fltarr(n_detx, n_dety) 
   dmask_fp2_exp=fltarr(n_detx, n_dety) 
   ; Kaya's code.

   for i=0, n_detx-1 do begin 
      for j=0, n_dety-1 do begin
         dmask_fp1(i,j) = (dmask(i,j) gt 0.0)? 1.0 : 0.0
         dmask_fp2(i,j) = (dmask(i+n_detx,j) gt 0.0)? 1.0 : 0.0
         dmask_fp1_exp(i,j) = (dmask(i,j) gt 0.0)? 0.0 : 1.0
         dmask_fp2_exp(i,j) = (dmask(i+n_detx,j) gt 0.0)? 0.0 : 1.0
      endfor 
   endfor

; here I modified original Kaya's code for percetage estimation 
; which gave 1/63^2 for zero SLP
index=where(dmask(0:63,*) gt 0.0,count)
fp1_pct=0.0
if(count gt 0) then fp1_pct=count*100.0/64.^2

index=where(dmask(64:127,*) gt 0.0,count)
fp2_pct=0.0
if(count gt 0) then fp2_pct=count*100.0/64.^2

fp1or2_pct=0.0
index=where (dmask_fp1 eq 0.0 or dmask_fp2 eq 0.0, count)
if(count gt 0) then fp1or2_pct=count*100.0/64.^2

efficiency=(100-fp1_pct)+(100-fp2_pct)


   if(silent ne 1) then begin
      !p.charsize=1.8
      print, 'FPA, FPB, 200-FPA-FPB, FPA.or.FPB [%]= ', fp1_pct, fp2_pct, (100-fp1_pct)+(100-fp2_pct), fp1or2_pct,format='(a,4f8.2)'
      contour, dmask_fp1, xpos_array, ypos_array, /cell_fill, xtit='DETX [mm]', ytit='DETY [mm]', $
               tit='FP1 '+'(PA = '+string(PAdeg, format='(f4.0)')+' deg)'
      oplot, [0], [0], psym=2    
      contour, dmask_fp2, xpos_array, ypos_array, /cell_fill, xtit='DETX [mm]', ytit='DETY [mm]', $
               tit='FP2 '+'(PA = '+string(PAdeg,format='(f4.0)')+' deg)'
      oplot, [0], [0], psym=2

      ;print,'obsid ',obsid[t]
      ;writefits, saveto+'map_fpA.obsid'+obsid[t]+'.'+key+'.fits',dmask_fp1_exp,hdr_float
      ;writefits, saveto+'map_fpB.obsid'+obsid[t]+'.'+key+'.fits',dmask_fp2_exp,hdr_float

   endif

   ; run over individual map for the best PA
   ; here we put individual observation GRID2 onto the GRID1 (Chandra map)
   for xx = 0, 63 do begin
      for yy = 0, 63 do begin
         ; take ra/dec of a given pixel on GRID2
         xy2ad, xx, yy, astr, ra, dec
         ; find the corresponding pixel on GRID1 (ra/dec => i/j)
         EULER, ra, dec, l, b, 1
         adxy, grid_hdr , l, b, i, j
         ; check limits
         if(i lt 0 or j lt 0 or i ge nwx or j ge nwy) then continue        
         map_cnt_pix(i,j)++
         map_fp1(i,j)+=dmask_fp1(xx,yy)
         map_fp2(i,j)+=dmask_fp2(xx,yy)
         map_fp1_exp(i,j)+=dmask_fp1_exp(xx,yy)
         map_fp2_exp(i,j)+=dmask_fp2_exp(xx,yy)
      endfor
   endfor

   if(0) then begin
   for xx = 0, nwx-1 do begin
      for yy = 0, nwy-1 do begin
         ; take ra/dec of a given pixel on GRID2
         xyad, grid_hdr, xx, yy, lon, lat
         ; find the corresponding pixel on GRID1 (ra/dec => i/j)
         EULER, lon, lat, ra, dec, 2
         ad2xy, ra, dec, astr, i, j
         if(i ge 0 and j ge 0 and i le 63 and j le 63) then begin
            fov_fpA=(dmask_fp1_exp(i,j) gt 0)? 1.0 : 0.0
            fov_fpB=(dmask_fp2_exp(i,j) gt 0)? 1.0 : 0.0
            map_cnt_fov_fpAB(xx,yy)+=(fov_fpA+fov_fpB)
            map_cnt_fov_fpA(xx,yy)+=fov_fpA
            map_cnt_fov_fpB(xx,yy)+=fov_fpB
         endif
      endfor
   endfor
   endif

endfor; end of observations cycle



; GRID3: Save grid for a given PA


map_fp1or2_pct=0.
map_fp1_pct=0.
map_fp2_pct=0.
map_total_pix=0.
map_total_open=0.
for i=0,nwx-1 do begin
   for j=0,nwy-1 do begin
      if(map_cnt_pix(i,j) eq 0) then continue
      map_total_pix+=map_cnt_pix(i,j)
      map_total_open++
      if(map_fp1(i,j) eq 0 or map_fp2(i,j) eq 0) then map_fp1or2_pct++
      if(map_fp1(i,j) gt 0) then map_fp1_pct+=map_fp1(i,j)
      if(map_fp2(i,j) gt 0) then map_fp2_pct+=map_fp2(i,j)
   endfor
endfor

fp1_pct=map_fp1_pct*100./map_total_pix
fp2_pct=map_fp2_pct*100./map_total_pix
map_eff=(100-fp1_pct)+(100-fp2_pct)
map_vis=map_fp1or2_pct*100./map_total_open


print, 'PA, FPA, FPB, 200-FPA-FPB, FPA.or.FPB [%]= ',PAdeg, fp1_pct, fp2_pct, map_eff, map_vis,format='(a,5f8.2)'
print,'total pix/open', map_total_pix,map_total_open


openw,1,saveto+'PA.'+key+'.dat',/append
printf,1, PAdeg, fp1_pct, fp2_pct, map_eff, map_vis,format='(a,5f8.2)'
close,1

; GRID1:
if(1) then begin
sxaddpar, grid_hdr,'PA', PAdeg
writefits, saveto+'map_fpA.'+key+'.fits',map_fp1,grid_hdr
writefits, saveto+'map_fpB.'+key+'.fits',map_fp2,grid_hdr
writefits, saveto+'map_exp_fpA.'+key+'.fits',map_fp1_exp,grid_hdr
writefits, saveto+'map_exp_fpB.'+key+'.fits',map_fp2_exp,grid_hdr
writefits, saveto+'map_fov_fpAB.'+key+'.fits',map_cnt_fov_fpAB,grid_hdr
writefits, saveto+'map_fov_fpA.'+key+'.fits',map_cnt_fov_fpA,grid_hdr
writefits, saveto+'map_fov_fpB.'+key+'.fits',map_cnt_fov_fpB,grid_hdr
writefits, saveto+'map_cnt_pix.'+key+'.fits',map_cnt_pix,grid_hdr
endif


endfor; end of PA cycle

;endfor; end of source's cycle

end
