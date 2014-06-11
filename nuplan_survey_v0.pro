pro nuplan_survey, pnt_file=pnt_file, map_center=map_center, slp_source_catalog_file=slp_source_catalog_file, silent=silent, key=key, pa_array=pa_array
Compile_Opt defint32

isGrid=0

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


; GRID3
; this is grid in galactic latitude, centered at Norma survey
; this grid is used only for high-level products like SLP coverage in %
gal=map_center
gwx=1.3d; 2 degrees wide in both directions
gwy=0.2d
hnwx=2; how many pixels in half of width
hnwy=1
nwx=hnwx*2+1; the full number of pixels (this was done to have center of the map in the central pixel) check maps carefully in ds9!
nwy=hnwy*2;+1
stepx=gwx/(nwx-1)
stepy=gwy/(nwy-1)
nn=nwx*nwy
grid_ra  = dblarr(nn); the position of the individual snapshot (ra)
grid_dec = dblarr(nn); the position of the individual snapshot (dec)
grid_x   = lonarr(nn)
grid_y   = lonarr(nn)
grid_fp1 = dblarr(nwx,nwy); this array is for SLP contamination in % (FP1)
grid_fp2 = dblarr(nwx,nwy); this array is for SLP contamination in % (FP2)
grid_eff_2d  = dblarr(nwx,nwy); PA for the every observation
grid_pa_2d  = dblarr(nwx,nwy); PA for the every observation
grid_pa_1d  = dblarr(nn); PA for the every observation
map_fp1=dblarr(nwx,nwy); the map for FP1
map_fp2=dblarr(nwx,nwy); the map for FP2
cnt=lonarr(nwx,nwy); exposure map

; make WCS coordinates for GRID3 (Note Cartesian projection -CAR)
  crpix=[fix(nwx/2.0d)+1,fix(nwy/2.0d)+1]
  cdelt=[-stepx,stepy]
  mkhdr, grid_hdr, 4, [nwx,nwy]
  make_astr,astr, DELTA = cdelt, CRPIX = crpix, CRVAL = [gal(0),gal(1)], CTYPE=['GLON-CAR','GLAT-CAR']
  putast, grid_hdr, astr, CD_TYPE=2

; GRID3

openw,1,'galactic_grid.reg'
for i=0,nwx-1 do begin
   for j=0,nwy-1 do begin
      ; make the corresponding lon/lat  
      l=gal(0)-gwx/2.0d + stepx*i
      b=gal(1)-gwy/2.0d + stepy*j
      ; convert lon/lat to ra/dec
      EULER, l, b, ra, dec, 2
      ; make running 1D index
      index=i*nwy+j
      ; survey individual pointings
      grid_ra(index)=ra
      grid_dec(index)=dec
      grid_x(index)=i
      grid_y(index)=j
      ;print,l,b,ra,dec
      printf,1,'galactic; circle(',l,', ',b,', 0.21)', FORMAT='(a,f12.4,a,f12.4,a)'
   endfor
endfor
close,1

; substitute grid with pointings file
if(~isGrid) then begin
   read_nustar_pointings,pnt_file,pnt_ra,pnt_dec,pnt_lon=pnt_lon,pnt_lat=pnt_lat
   n_pnt=n_elements(pnt_ra)
endif else begin
   n_pnt=nn
endelse

koeff=0.0d; not used
!p.multi=[0,2,1]


; read catalog of sources
src_file=root+slp_source_catalog_file
read_source_catalog,src_file,src_name,src_ra,src_dec,src_flag,index=index

; save source list used in this run
openw,1,saveto+'sources.'+key+'.dat'
for ii=0, n_elements(src_name)-1 do printf,1,src_name(ii),',',src_ra(ii),', ',src_dec(ii),',',src_flag(ii)
close,1

;for l=0,n_elements(src_name)-1 do begin
;key=STRCOMPRESS(src_name[l],/REMOVE_ALL)
;RA_SRC_IN=[src_ra[l]]
;DEC_SRC_IN=[src_dec[l]]

RA_SRC_IN=src_ra(index)
DEC_SRC_IN=src_dec(index)

RA_SRC = RA_SRC_IN*Deg2Rad
DEC_SRC = DEC_SRC_IN*Deg2Rad
Flux = replicate(1., n_elements(RA_SRC))
n_src = n_elements(RA_SRC) 

file_delete,saveto+'PA.'+key+'.dat',/ALLOW_NONEXISTENT
; set PA array:

for r=0,n_elements(pa_array)-1 do begin
   PAdeg=-180.0-pa_array[r]
   PA=PAdeg*Deg2Rad             ; in radians

   cnt(*,*)=0
   map_fp1(*,*)=0.
   map_fp2(*,*)=0.

; running over survey individual pointings
for t=0,n_pnt-1 do begin
   ; take fixed PA angle
   ;PAdeg=2.0

   if(isGrid) then begin
      this_pnt_ra=grid_ra(t)
      this_pnt_dec=grid_dec(t)
   endif else begin
      this_pnt_ra=pnt_ra(t)
      this_pnt_dec=pnt_dec(t)
   endelse
   
   ; make header for individual observation
   make_astr_nustar, this_pnt_ra, this_pnt_dec, PAdeg, astr=astr

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
   ; Kaya's code.
   for i=0, n_detx-1 do begin 
      for j=0, n_dety-1 do begin 
         if dmask(i,j) gt 0.0 then begin 
            dmask_fp1(i,j) = dmask(i, j)/dmask(i,j)
         endif else begin 
            dmask_fp1(i,j) = 0.0 
         endelse 
         if dmask(i+n_detx,j) gt 0.0 then begin
            dmask_fp2(i,j) = dmask(i+n_detx, j)/dmask(i+n_detx,j)
         endif else begin
            dmask_fp2(i,j) = 0.0
         endelse
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

if(isGrid) then begin
   grid_pa_1d(t)=PAdeg
   grid_eff_2d(grid_x(t),grid_y(t))=efficiency
   grid_pa_2d(grid_x(t),grid_y(t))=PAdeg
   grid_fp1(grid_x(t),grid_y(t))=fp1_pct
   grid_fp2(grid_x(t),grid_y(t))=fp2_pct
endif else begin
endelse


   if(silent ne 1) then begin
      !p.charsize=1.8
      print, 'FPA, FPB, 200-FPA-FPB, FPA.or.FPB [%]= ', fp1_pct, fp2_pct, (100-fp1_pct)+(100-fp2_pct), fp1or2_pct,format='(a,4f8.2)'
      contour, dmask_fp1, xpos_array, ypos_array, /cell_fill, xtit='DETX [mm]', ytit='DETY [mm]', $
               tit='FP1 '+'(PA = '+string(PAdeg, format='(f4.0)')+' deg)'
      oplot, [3], [3], psym=2    
      contour, dmask_fp2, xpos_array, ypos_array, /cell_fill, xtit='DETX [mm]', ytit='DETY [mm]', $
               tit='FP2 '+'(PA = '+string(PAdeg,format='(f4.0)')+' deg)'
      oplot, [3], [3], psym=2
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
         ; "exposure" map of GRID1
        
         cnt(i,j)++
         if(xx eq 38 and yy eq 38) then cnt(i,j)+=100
         map_fp1(i,j)+=dmask_fp1(xx,yy)
         map_fp2(i,j)+=dmask_fp2(xx,yy)
      endfor
   endfor

endfor; end of observations cycle


; GRID3: Save grid for a given PA

if(isGrid) then begin
;sxaddpar, grid_hdr,'PA', PAdeg
;writefits,'grid_fp1.'+key+'.fits',grid_fp1,grid_hdr
;writefits,'grid_fp2.'+key+'.fits',grid_fp2,grid_hdr
;writefits,'grid_pa.'+key+'.fits',grid_pa_2d,grid_hdr
;writefits,'grid_eff.'+key+'.fits',grid_eff_2d,grid_hdr
;writefits, 'cnt.'+key+'.fits',cnt,h
;save,filename='grid.'+key+'.sav',grid_lon,grid_lat,grid_ra,grid_dec,grid_pa_1d,grid_pa_2d,grid_eff_2d,grid_hdr,PAdeg
endif else begin
endelse

map_fp1or2_pct=0.
map_fp1_pct=0.
map_fp2_pct=0.
map_total_pix=0.
map_total_open=0.
for i=0,gwx-1 do begin
   for j=0,gwy-1 do begin
      if(cnt(i,j) eq 0) then continue
      map_total_pix+=cnt(i,j)
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

openw,1,saveto+'PA.'+key+'.dat',/append
printf,1, PAdeg, fp1_pct, fp2_pct, map_eff, map_vis,format='(a,5f8.2)'
close,1

; GRID1:
if(1) then begin
sxaddpar, grid_hdr,'PA', PAdeg
writefits, 'map_fp1.'+key+'.fits',map_fp1,grid_hdr
writefits, 'map_fp2.'+key+'.fits',map_fp2,grid_hdr
writefits, 'cnt.'+key+'.fits',cnt,grid_hdr
endif

endfor; end of PA cycle

;endfor; end of source's cycle

end
