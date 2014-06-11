;;
;; Makes mosaic resampling pixel counts over the new grid
;;

; Example:
;nustar_mosaic, 266.45958, -28.821944, cdelt=3, width=800, key='03to10keV', $
;               root='/Data/NuSTAR/data/mini-survey', $
;               obsid=['40010003001A', '40010005001A'], $
;               evtpath='event_cl_new_slp_v1'


;; 
;;
;;

pro mosaic, center_ra, center_dec, width=width, rota=rota, nn=nn, e1=e1, e2=e2, $
            key=key, cdelt=cdelt, root=root, obsid=obsid, evtpath=evtpath, $
            devmax=devmax, expocut=expocut, psf=psf, psf_ra=psf_ra, psf_dec=psf_dec, offhisto_path=offhisto_path

dr = !PI/180.

if(n_elements(psf) eq 0) then psf=0
if(n_elements(offhisto_path) eq 0) then offhisto_path='./products'
if(n_elements(psf_ra) eq 0) then psf_ra=center_ra
if(n_elements(psf_dec) eq 0) then psf_dec=center_dec

if(n_elements(devmax) eq 0) then devmax=-1.0
if(n_elements(cdelt) eq 0) then cdelt=2.0
if(n_elements(rota) eq 0) then rota=0.0
if(n_elements(nn) eq 0) then nn=3
if(n_elements(e1) eq 0) then e1=3.0
if(n_elements(e2) eq 0) then e2=10.0
if(n_elements(evtpath) eq 0) then evtpath='event_cl'
if(n_elements(expocut) eq 0) then expocut=10000.0

; width of the mosaic in seconds
if(n_elements(width) eq 0) then width=800.0


if(n_elements(root) eq 0) then begin
   print,'Give me root directory! root=?'
   return
endif

if(n_elements(root) eq 0) then begin
   print,"Give me list of observations! e.g. obsid=['40010003001A', '40010005001A', '40010006001A', '40010005001B', '40010006001B']"
   return
endif

keyout=key


w=Fix(width/cdelt)+1
s2exp=DBLARR(w,w)
s2cnt=DBLARR(w,w)
s2d=DBLARR(w,w)

crpix=[w/2.,w/2.]
crval=[center_ra, center_dec]
rot_mat = [ [ cos(rota*dr), -sin(rota*dr)], $ ;Rotation matrix, PA doesn't matter
            [ sin(rota*dr),  cos(rota*dr)] ] 
make_astr,astr_map, CD=rot_mat, DELT = [cdelt/3600., cdelt/3600.], CRPIX = crpix, CRVAL = crval, $
           RADECSYS = 'FK5', EQUINOX = 2000.0
mkhdr, hdr_map, 4, [w,w] 
putast, hdr_map, astr_map, CD_TYPE=2

nn2=nn^2

s2exp[*,*]=0.0
s2cnt[*,*]=0.0
s2d[*,*]=0.0

for di=0L,n_elements(obsid)-1 do begin
   obs=obsid[di]

;print,'Make sky...'
obsid_s = STRMID( obs, 0,11)
module = STRMID( obs, 11, 1)
if(module eq 'A') then imod=0
if(module eq 'B') then imod=1


make_sky,root=root, obsid=obsid_s, silent=silent, $
         evtpath=evtpath, e1=e1, e2=e2, key=key,  sky=sky, hdr=hdr_sky, delta=delta, PA=PA

s1d=Double(sky[*,*,imod])
sz=size(sky)
sx=sz[1]
sy=sz[2]

;; load PSF
   off_histo_name=offhisto_path+'/nu'+obs+'01_offaxishisto.fits'
   if(psf ne 0) then begin
   if(file_test(off_histo_name)) then begin
      print, 'Reading ',off_histo_name
      adxy,hdr_sky,psf_ra,psf_dec,psf_x,psf_y
      print, 'PA: ',PA
      print,'psf_x, psf_y: ',psf_x, psf_y
      CreateEffPSF, psf_x, psf_y, PA, psf_x, psf_y, off_histo_name, PSFeff=PSFeff, hdr_psf=hdr_psf, $
                    astr_psf=astr_psf, hdr_sky=hdr_sky, arcminPerPixPSF=arcminPerPixPSF, psf_koeff=1.0
      ;print,'arcminPerPixPSF:',arcminPerPixPSF
      sz=size(PSFeff)
      psf_sx=sz[1]
      psf_sy=sz[2]
      s1psf = PSFeff/total(PSFeff) ;*obsid_live[di]/live_total ; normalize the psf to total=1, and weight according to exposure
      writefits,obs+'.s1psf.fits',s1psf,hdr_psf
   endif else begin
      ;print, '*** Warning *** Does not exist: ',off_histo_name,' *** skip ***'
   endelse
   endif
;; end of loading PSF



expo_file=root+'/'+obsid_s+'/'+evtpath+'/nu'+obs+'01_ex.img'
if(~file_test(expo_file)) then begin
   print,'Exposure is not found:'
   print, expo_file
   print, 'run $NUPLAN/nuselect.pl in '+root+'/'+obsid_s+'/'+evtpath
   continue
endif
s1exp=readfits(expo_file,hdr_expo,/SILENT)
index=where(s1exp lt expocut)
s1exp(index)=-1e9

; Accumulate sky
for xi=0L,sx-1 do for yi=0L,sy-1 do begin
   xyad,hdr_sky,xi,yi,ra,dec
   if (devmax gt 0.0) then begin
   dist=sphdist(center_ra, center_dec, ra, dec, /DEGREES)*3600.
   if(dist gt devmax) then continue
   endif

   for ii=1L,nn do for jj=1L,nn do begin
      x=Double(xi)-0.5+(Double(ii)-0.5)/Double(nn)
      y=Double(yi)-0.5+(Double(jj)-0.5)/Double(nn)
      xyad,hdr_sky,x,y,ra,dec
      adxy,hdr_map,ra,dec,mapx,mapy
      mapi=Round(mapx,/L64)
      mapj=Round(mapy,/L64)

      ; Add fraction of the pixel to the base sky
      if (mapi ge 0 and mapi lt w and mapj ge 0 and mapj lt w) then begin
         if(s1exp[xi,yi] gt 0.0) then begin
            s2exp[mapi,mapj]+=s1exp[xi,yi]/nn2 
            s2d[mapi,mapj]+=(s1d[xi,yi])/nn2
         endif
      endif
   endfor
endfor


endfor ; obsid cycle


s2d_out=s2d
writefits,'map.'+keyout+'.cts.fits',s2d,hdr_map
; Convert summed sky back to normal units
for xi=0L,w-1 do for yi=0L,w-1 do begin
   if(s2exp[xi,yi] gt 0.0) then begin
      s2d_out[xi,yi]=s2d[xi,yi]/s2exp[xi,yi]
   endif else begin
      s2d_out[xi,yi]=0
   endelse
endfor
writefits,'map.'+keyout+'.exp.fits',s2exp,hdr_map
writefits,'map.'+keyout+'.flx.fits',s2d_out,hdr_map


end
