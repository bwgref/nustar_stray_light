;;
;; Makes mosaic resampling pixel counts over the new grid
;; Needs exposure map after nuselect.pl
;; step 1:
;;        edit nuselect.pl and specify your energy bands
;;        go to your 'event_cl' and run nuselect.pl (it makes e-bands
;;        dir with exposure maps and others)
;; step 2:
;;        run this code keeping the same key as nuselect.pl uses
;;

; Example:
;pro runme
;nustar_mosaic, 266.45958, -28.821944, cdelt=cdelt, width=800, key='03to10keV', $
;               root='/Data/NuSTAR/data/mini-survey', $
;               obsid=['40010003001A', '40010005001A'], $
;               evtpath='event_cl_new_slp_v1'
;end

;; 
;;
;;

function annulus_mean_flux,ra,dec,Rmin,Rmax,map,hdr,expo=expo,$
                           sx=sx,sy=sy,bx=bx,by=by,smax=smax,bmax=bmax,delta=delta
  
  npix=Ceil(3.14159*(Rmax/(3600*abs(delta[1]))+5)^2) ; --> +5 pixels
  sx=LONARR(npix)
  sy=LONARR(npix)
  bx=LONARR(npix)
  by=LONARR(npix)
  smax=0L
  bmax=0L

  sz=size(map)
  wx=sz[1]
  wy=sz[2]
  ;bkg=LONARR(sx,sy)
  totpix=0L
  totcnt=0.0
  for xi=0L,wx-1 do for yi=0L,wy-1 do begin
     xyad,hdr,xi,yi,xi_ra,xi_dec
     if not (expo[xi,yi] gt 0) then continue
     dist=sphdist(ra, dec, xi_ra, xi_dec, /DEGREES)*3600.
     if(dist gt Rmax) then continue


     if(smax ge npix or bmax ge npix) then begin
        print,'abort',smax,bmax,npix
        stop
     endif
     if(dist lt Rmin) then begin
        sx[smax]=xi
        sy[smax]=yi
        smax++
     endif else begin
        bx[bmax]=xi
        by[bmax]=yi
        bmax++
     endelse

     if keyword_set(expo) then if(expo[xi,yi] eq 0) then continue
     ;bkg[xi,yi]=map[xi,yi]
     totcnt+=map[xi,yi]/expo[xi,yi]
     totpix++
  endfor
  ;writefits,'bkg.fits',bkg,hdr

  return, Double(totcnt)/Double(totpix)
end


pro stack_catalog_obsid, catalog=catalog, width=width, rota=rota, nn=nn, e1=e1, e2=e2, $
            key=key, cdelt=cdelt, root=root, obsid=obsid, evtpath=evtpath, $
            devmax=devmax, expocut=expocut, out=out, Rmin=Rmin, Rmax=Rmax

center_ra=0.0
center_dec=0.0
dr = !PI/180.

if(n_elements(psf) eq 0) then psf=0
if(n_elements(offhisto_path) eq 0) then offhisto_path='./products'

if(n_elements(devmax) eq 0) then devmax=600.0
if(n_elements(cdelt) eq 0) then cdelt=2.0
if(n_elements(rota) eq 0) then rota=0.0
if(n_elements(nn) eq 0) then nn=3
if(n_elements(e1) eq 0) then e1=3.0
if(n_elements(e2) eq 0) then e2=10.0
if(n_elements(evtpath) eq 0) then evtpath='event_cl'
if(n_elements(out) eq 0) then out=''
if(n_elements(expocut) eq 0) then expocut=10000.0

; width of the mosaic in seconds
if(n_elements(width) eq 0) then width=800.0


if(n_elements(root) eq 0) then begin
   print,'Give me root directory! root=?'
   return
endif

if(n_elements(root) eq 0) then begin
   print,"Give me obsid! e.g. obsid='40010003001A'"
   return
endif

keyout=key+'.'+obsid
if(out ne '') then keyout=key+'.'+obsid+'.'+out


w=Fix(width/cdelt)+1
s2psf=DBLARR(w,w)
s2exp=DBLARR(w,w)
s2cnt=DBLARR(w,w)
s2d=DBLARR(w,w)
s2d2=DBLARR(w,w)
s2e=DBLARR(w,w)

crpix=[w/2.,w/2.]
crval=[0.0, 0.0]
rot_mat = [ [ cos(rota*dr), -sin(rota*dr)], $ ;Rotation matrix, PA doesn't matter
            [ sin(rota*dr),  cos(rota*dr)] ] 
make_astr,astr_map, CD=rot_mat, DELT = [cdelt/3600., cdelt/3600.], $
          CRPIX = crpix, CRVAL = crval, $
          RADECSYS = 'FK5', EQUINOX = 2000.0
mkhdr, hdr_map, 4, [w,w] 
putast, hdr_map, astr_map, CD_TYPE=2

nn2=nn^2

s2psf[*,*]=0.0
s2exp[*,*]=0.0
s2cnt[*,*]=0.0
s2d[*,*]=0.0
s2d2[*,*]=0.0
s2e[*,*]=0.0

read_nxy_catalog,catalog,src_name,src_ra,src_dec


;for di=0L,n_elements(obsid)-1 do begin
obs=obsid
print,'ObsID ',obs

obsid_s = STRMID( obs, 0,11)
module = STRMID( obs, 11, 1)
if(module eq 'A') then imod=0
if(module eq 'B') then imod=1

make_sky,root=root, obsid=obsid_s, silent=silent, evtpath=evtpath,$
         e1=e1, e2=e2, key=key,  sky=sky, hdr=hdr_sky,  delta=delta, PA=PA, $
         center_ra=center_ra, center_dec=center_dec

s1d=Double(sky[*,*,imod])
sz=size(sky)
dx=sz[1]
dy=sz[2]


expo_file=root+'/'+obsid_s+'/'+evtpath+'/nu'+obs+'01_ex.img'
if(~file_test(expo_file)) then begin
   print,'Exposure not found:'
   print, expo_file
   return
endif
print,'Read: ',expo_file

s1exp=readfits(expo_file,hdr_expo,/SILENT)
sz=size(s1exp)
sxe=sz[1]
sye=sz[2]
index=where(s1exp lt expocut)
s1exp(index)=0.0

Rmin=30.
Rmax=100.
FOVmax=7.0 ; arcmin
livetime_total=0.0
for ss=0L,N_ELEMENTS(src_ra)-1 do begin
   ;print,'1: src ',ss,' / ',N_ELEMENTS(src_ra),' ',obs
   ra=src_ra[ss]
   dec=src_dec[ss]
   offset=sphdist(center_ra, center_dec, ra, dec, /DEGREES)*60.
   if(offset gt FOVmax) then continue

   ;; Read exposure
   adxy,hdr_expo,src_ra[ss],src_dec[ss],xi,yi
   if(xi lt 0 or yi lt 0 or xi ge sxe or yi ge sye) then begin
      print,'Wrong exposure position'
      stop
   endif
   livetime=s1exp[Round(xi,/L64),Round(yi,/L64)]
   print,livetime
   if not (livetime gt 0.0) then continue
   ;; Load PSF
   off_histo_name='./products/'+obsid+'_'+src_name[ss]+'_offaxishisto.fits'
   if not (file_test(off_histo_name)) then continue
   print,'Reading ',off_histo_name
   adxy,hdr_sky,src_ra[ss],src_dec[ss],src_x,src_y
   spawn, "ls "+off_histo_name
   CreateEffPSF, src_x, src_y, PA, src_x, src_y, off_histo_name, PSFeff=PSFeff, hdr_psf=hdr_psf, $
                 astr_psf=astr_psf, hdr_sky=hdr_sky, arcminPerPixPSF=arcminPerPixPSF, psf_koeff=1.0
   sz=size(PSFeff)
   psf_sx=sz[1]
   psf_sy=sz[2]
   s1psf = PSFeff/total(PSFeff)*livetime
   livetime_total+=livetime
   writefits,obsid+'_'+src_name[ss]+'.s1psf.fits',s1psf,hdr_psf
   ;; end of PSF loading part

   make_astr,astr_src, CD=rot_mat, DELT = [cdelt/3600., cdelt/3600.], $
             CRPIX = crpix, CRVAL = [ra,dec], RADECSYS = 'FK5', EQUINOX = 2000.0
   mkhdr, hdr_src, 4, [w,w] 
   putast, hdr_src, astr_src, CD_TYPE=2

   adxy,hdr_sky,ra,dec,x,y
   xi=Round(x,/L64)
   yi=Round(y,/L64)
   if not (xi ge 0 and xi lt dx and yi ge 0 and yi lt dy) then continue
   if not (s1exp[xi,yi] gt 0) then continue
   bkg=annulus_mean_flux(ra,dec,Rmin,Rmax,sky,hdr_sky,expo=s1exp,delta=delta,sx=sx,sy=sy,bx=bx,by=by,$
                         smax=smax,bmax=bmax)
   if(bmax*smax eq 0) then continue
   ;print,'2: src ',ss,' / ',N_ELEMENTS(src_ra)

   brate=0.0
   btot=0.0
   for xx=0L,bmax do begin
      btot+=s1exp[bx[xx],by[xx]]
      brate+=s1d[bx[xx],by[xx]]
   endfor
 
   for xx=0L,bmax do begin
      xi=bx[xx]
      yi=by[xx]

      for ii=1L,nn do for jj=1L,nn do begin
         x=Double(xi)-0.5+(Double(ii)-0.5)/Double(nn)
         y=Double(yi)-0.5+(Double(jj)-0.5)/Double(nn)

         xyad,hdr_sky,x,y,xi_ra,xi_dec
         adxy,hdr_src,xi_ra,xi_dec,m,n

         mi=Round(m,/L64)
         ni=Round(n,/L64)

         if (mi ge 0 and mi lt w and ni ge 0 and ni lt w) then begin
               s2exp[mi,ni]+=s1exp[xi,yi]/nn2 
               s2d[mi,ni]+=(s1d[xi,yi])/nn2
               s2d2[mi,ni]+=(s1d[xi,yi]-brate*s1exp[xi,yi])/nn2
               s2cnt[mi,ni]++
         endif
      endfor
   endfor
   
   ;print,'3: src ',ss,' / ',N_ELEMENTS(src_ra)

   for xx=0L,smax do begin
      xi=sx[xx]
      yi=sy[xx]
  
      for ii=1L,nn do for jj=1L,nn do begin
         x=Double(xi)-0.5+(Double(ii)-0.5)/Double(nn)
         y=Double(yi)-0.5+(Double(jj)-0.5)/Double(nn)

         xyad,hdr_sky,x,y,xi_ra,xi_dec
         adxy,hdr_src,xi_ra,xi_dec,m,n

         mi=Round(m,/L64)
         ni=Round(n,/L64)

         if (mi ge 0 and mi lt w and ni ge 0 and ni lt w) then begin
               s2exp[mi,ni]+=s1exp[xi,yi]/nn2
               s2d[mi,ni]+=(s1d[xi,yi])/nn2
               s2d2[mi,ni]+=(s1d[xi,yi]-brate*s1exp[xi,yi])/nn2
               s2e[mi,ni]+=(brate*s1exp[xi,yi])/nn2
               s2cnt[mi,ni]++
         endif
      endfor
   endfor

;; project PSF
   for xi=0L,psf_sx-1 do for yi=0L,psf_sy-1 do begin
      for ii=1L,nn do for jj=1L,nn do begin
         x=Double(xi)-0.5+(Double(ii)-0.5)/Double(nn)
         y=Double(yi)-0.5+(Double(jj)-0.5)/Double(nn)
         
         xyad,hdr_psf,x,y,xi_ra,xi_dec
         adxy,hdr_src,xi_ra,xi_dec,m,n
      
         mi=Round(m,/L64)
         ni=Round(n,/L64)
         
         if (mi ge 0 and mi lt w and ni ge 0 and ni lt w) then begin
            s2psf[mi,ni]+=s1psf[xi,yi]/nn2
         endif
      endfor   
   endfor

;if (ss gt 0) then continue

endfor




;endfor ; obsid cycle



s2d_out=s2d
s2d2_out=s2d2
s2e_out=s2e
writefits,'map.'+keyout+'.cts.fits',s2d,hdr_map
print,'Convert summed sky back to normal units'
for xi=0L,w-1 do for yi=0L,w-1 do begin
   if(s2exp[xi,yi] gt 0.0) then begin
      s2d_out[xi,yi]=s2d[xi,yi]/s2exp[xi,yi]
      s2d2_out[xi,yi]=s2d2[xi,yi]/s2exp[xi,yi]
      s2e_out[xi,yi]=sqrt(s2e[xi,yi])/s2exp[xi,yi]
   endif else begin
      s2d_out[xi,yi]=0
      s2d2_out[xi,yi]=0
      s2e_out[xi,yi]=0
   endelse
endfor
writefits,'map.'+keyout+'.exp.fits',s2exp,hdr_map
writefits,'map.'+keyout+'.err.fits',s2e_out,hdr_map
writefits,'map.'+keyout+'.flx.fits',s2d_out,hdr_map
writefits,'map.'+keyout+'.flx-resid.fits',s2d2_out,hdr_map

writefits,'map.'+keyout+'.psf.fits',s2psf/livetime_total,hdr_map

end
