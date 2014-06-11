; OAA [radian] 
; AZIMUTH: Azimuthal angle [radian]                                                                                                              
; PA: pointing position angle [radian] which describes the rotation of the
; DETECTOR coordinates wrt SKY

FUNCTION LEAKAGE_MAP, OAA, AZIMUTH, PA 

; Set up 

my_azimuth = !PI/2. - AZIMUTH ; Peter's (Mao?) original line 
;my_azimuth = AZIMUTH + !PI  ; equation that works for 1E1740 and Sgr A* 

DET2FB = complex(cos(2.*!PI/3.), sin(2.*!PI/3.))
SKY2DET = complex(cos(PA), sin(PA)) ; Original formula from Peter's code 
;SKY2DET = complex(cos(!PI-PA), sin(!PI-PA)) ; modified to match
;(returned back by Roman)
SKY2FB = SKY2DET * DET2FB 

cbox_real = [1, 1, -1, -1, 1] 
cbox_im = [1, -1, -1, 1, 1] 
cbox = complex(cbox_real, cbox_im) 

fpdet = cbox * 20 * DET2FB

det0 = complex(cbox*9.5+9.75, 9.75) * DET2FB 
det3 = complex(cbox*9.5+9.75, -9.75) * DET2FB 

iFL = 1./10.015 ; [1/meter] inverse focal length (converts mm to mrad) 
hgap = 0.15 ; detector half gap [mm] 

; Bench 

bench_z = 10.2e3 
bench_r_real = [-64,-176,-237,-234,-118,173,334,626,747,748,689,570,-64]
bench_r_im = [-240,-140,-44,92,205,588,588,205,92,-44,-140,-240,-240]
;bench_r_real = [-64., -176., -237., -234., 173., 334., 747., 748., 689., 570., -64.] 
;bench_r_im = [-240., -140., -44., 92., 588., 588., 92., -44., -140., -240., -240.] 
bench_r = complex(bench_r_real, bench_r_im)
;bench_r = complex(bench_r_real, -bench_r_im) ; Peter's suggested sign flip to reproduce 1E1740 pattern (7/30)

bench_shadow = bench_r + (OAA*bench_z) * complex(cos(my_azimuth), sin(my_azimuth)) * SKY2FB 
bench_th = atan(bench_r) 
bench_phi = atan(abs(bench_r)/bench_z) 

; Aperture 1 

ap1_z = 833.2 ; 
ap1_r_array = findgen(101)*(2.*!PI/100.)+0.
ap1_r = (58./2.) * complex(cos(ap1_r_array), sin(ap1_r_array)) 
ap1_shadow = ap1_r + (OAA * ap1_z) * complex(cos(my_azimuth), sin(my_azimuth)) * SKY2FB
ap1_th = atan(ap1_r) 
ap1_phi = atan(abs(ap1_r)/ap1_z) 

; bitmask
d0mask = fltarr(64,64)
d1mask = fltarr(64,64)
xpos_array = fltarr(64)
ypos_array = fltarr(64)

ap1_x = real_part(ap1_shadow/DET2FB)
ap1_y= imaginary(ap1_shadow/DET2FB)
bench1_x =  real_part(bench_shadow/DET2FB)
bench1_y = imaginary(bench_shadow/DET2FB)
bench2_x =  real_part((bench_shadow-508.)/DET2FB)
bench2_y =imaginary((bench_shadow-508.)/DET2FB)


;; time = systime(1)
;; for xx = 0, 63 do begin 
;;    for yy = 0, 63 do begin 
;; ;      q_x= sign(xx+1 - 32.5)                                           ;
;; ;      q_y = sign(yy+1 - 32.5) 
      
;;       xpos = (xx+1-32.5)*0.6048 + hgap*sign(xx+1-32.5)
;;       ypos = (yy+1-32.5)*0.6048 + hgap*sign(yy+1-32.5)

;;       d0mask(xx,yy) = InsidePolygon(xpos, ypos, ap1_x, ap1_y) * (1.-InsidePolygon(xpos, ypos, bench1_x, bench1_y))
;;       d1mask(xx,yy) = InsidePolygon(xpos, ypos, ap1_x, ap1_y) * (1.-InsidePolygon(xpos, ypos, bench2_x, bench2_y)) 

;; ;      xpos_array(xx) = xpos                          
;; ;      ypos_array(yy) = ypos 
;;    endfor
;; endfor
;; print, systime(1) - time

;time = systime(1)
xx = findgen(64*64) mod 64
yy = floor(findgen(64*64) / 64.)
xpos = (xx+1-32.5)*0.6048 + hgap*sign(xx+1-32.5)
ypos = (yy+1-32.5)*0.6048 + hgap*sign(yy+1-32.5)
d0mask =  insidepolygon_vect(xpos, ypos, ap1_x, ap1_y) * (1.-InsidePolygon_vect(xpos, ypos, bench1_x, bench1_y) )
d0mask = reform(d0mask, [64, 64])

d1mask = insidepolygon_vect(xpos, ypos, ap1_x, ap1_y) * (1.-InsidePolygon_vect(xpos, ypos, bench2_x, bench2_y)) 
d1mask = reform(d1mask, [64, 64])
;print, systime(1) - time

;stop




;!p.multi=[0,2,1] 

;contour, d0mask, xpos_array, ypos_array, /cell_fill, tit='FP0 mask in DET coordinates', xtit='DET x [mm]', ytit='DET y [mm]'

;contour, d1mask, xpos_array, ypos_array, /cell_fill, tit='FP1 mask in DET coordinates', xtit='DET x [mm]', ytit='DET y [mm]'

return, [d0mask, d1mask] 


; Plots 

!p.multi=[0, 2, 3]

; Sky coordinates 

plot, real_part(-bench_shadow/SKY2FB*iFL), imaginary(bench_shadow/SKY2FB*iFL), tit='FP0 in SKY coordinates', xtit='SKY RA [mm]', ytit='SKY DEC [mm]', xr=[-10, 10], yr=[-10,10]

oplot, real_part(fpdet/SKY2FB*iFL), imaginary(fpdet/SKY2FB*iFL)

oplot, real_part(det0/SKY2FB*iFL), imaginary(det0/SKY2FB*iFL)

oplot, real_part(det3/SKY2FB*iFL), imaginary(det3/SKY2FB*iFL)

oplot, real_part(ap1_shadow/SKY2FB*iFL), imaginary(ap1_shadow/SKY2FB*iFL)

plot, real_part((bench_shadow-508.)/SKY2FB*iFL), imaginary((bench_shadow-508.)/SKY2FB*iFL), tit='FP1 in SKY coordinates', xtit='SKY RA [mm]', ytit='SKT DEC [mm]'

oplot, real_part(fpdet/SKY2FB*iFL), imaginary(fpdet/SKY2FB*iFL)

oplot, real_part(det0/SKY2FB*iFL), imaginary(det0/SKY2FB*iFL)

oplot, real_part(det3/SKY2FB*iFL), imaginary(det3/SKY2FB*iFL)

oplot, real_part(ap1_shadow/SKY2FB*iFL), imaginary(ap1_shadow/SKY2FB*iFL)

; Detector coordinates 

plot, real_part(bench_shadow/SKY2FB*iFL), imaginary(bench_shadow/DET2FB), tit='FP0 in DET coordinates', xtit='DET RA [mm]', ytit='DET DEC [mm]'

oplot, real_part(fpdet/SKY2FB*iFL), imaginary(fpdet/DET2FB)

oplot, real_part(det0/SKY2FB*iFL), imaginary(det0/DET2FB)

oplot, real_part(det3/SKY2FB*iFL), imaginary(det3/DET2FB)

oplot, real_part(ap1_shadow/SKY2FB*iFL), imaginary(ap1_shadow/DET2FB)

plot, real_part((bench_shadow-508.)/SKY2FB*iFL), imaginary((bench_shadow-508.)/DET2FB), tit='FP1 in DET coordinates', xtit='DET RA [mm]', ytit='DET DEC [mm]'

oplot, real_part(fpdet/SKY2FB*iFL), imaginary(fpdet/DET2FB)

oplot, real_part(det0/SKY2FB*iFL), imaginary(det0/DET2FB)

oplot, real_part(det3/SKY2FB*iFL), imaginary(det3/DET2FB)

oplot, real_part(ap1_shadow/SKY2FB*iFL), imaginary(ap1_shadow/DET2FB)

; FB coordinate plot 

plot, real_part(bench_shadow), imaginary(bench_shadow), tit='FP0 in FB coordinates', xtit='FB RA [mm]', ytit='FB DEC [mm]'

oplot, real_part(fpdet), imaginary(fpdet)

oplot, real_part(det0), imaginary(det0)

oplot, real_part(det3), imaginary(det3)

oplot, real_part(ap1_shadow), imaginary(ap1_shadow)

plot, real_part(bench_shadow), imaginary(bench_shadow), tit='FP1 in FB coordinates', xtit='FB RA [mm]', ytit='FB DEC [mm]'

oplot, real_part(fpdet+508.), imaginary(fpdet+508.)

oplot, real_part(det0+508.), imaginary(det0+508.)

oplot, real_part(det3+508.), imaginary(det3+508.)

oplot, real_part(ap1_shadow+508.), imaginary(ap1_shadow+508.)

return, [d0mask, d1mask] 


END 
