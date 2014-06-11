; NUPLAN_INSTRMAP adopted from Dan's nuskybgd code
;
; Creates an instrument map excluding RAW pixels contained in bad pixel files.
; The 'files' argument can be ignored if only the CALDB bad pixel list is used.
; Otherwise, files gives the name (and full path) to additional bad pixel lists --
; can be a single string or array of strings.


pro nuplan_instrmap,ab,clobber=clobber,pixmapx=pixmapx,pixmapy=pixmapy,detnum=detnum

if not keyword_set(clobber) then clobber=0

;grade weighting from NGC253 002 obs.
gw=[1.00000,    0.124902,    0.117130,    0.114720,    0.118038,   0.0114296,$
    0.0101738,   0.0113617,   0.0122017,   0.0157910,   0.0144079,   0.0145691,$
    0.0149934,  0.00165462,  0.00194312,  0.00156128,  0.00143400,  0.00210433,$
   0.00180735,  0.00140006,  0.00169704,  0.00189220,  0.00160371,  0.00150188,$
   0.00168007, 0.000296983, 0.000364864]

if(getenv('CALDB') eq '') then begin
   print,'*** ERROR ***'
   print,'Hey, where is NuSTAR $CALDB variable?'
   stop
endif

caldb=getenv('CALDB')+'/'
nuauxil=getenv('NUPLAN_AUXIL')
if (nuauxil eq '') then begin
   print,'*** ERROR ***'
   print,'1) make temporal working directory somewhere'
   print,'2) Set up NUPLAN_AUXIL variable in your shell'
   stop
endif
auxildir=nuauxil+'/'

dirinstrmap=caldb+'data/nustar/fpm/bcf/instrmap/'
dirpixpos=caldb+'data/nustar/fpm/bcf/pixpos/'
dirbadpix=caldb+'data/nustar/fpm/bcf/badpix/'

fits_read,dirinstrmap+'nu'+ab+'instrmap20100101v003.fits',instrmap,header

; check if images already exist
if not file_test(auxildir+'pixmapx'+ab+'.fits') or clobber then begin
    file=dirpixpos+'nu'+ab+'pixpos20100101v005.fits'
    pixmapx=lonarr(360,360)
    pixmapy=lonarr(360,360)
    pixmapx[*,*]=-1
    pixmapy[*,*]=-1
    detnum=lonarr(360,360)
    detnum[*,*]=-1
    allpdf=fltarr(360,360)
    for idet=0,3 do begin
      pixpos=mrdfits(file,idet+1,hh,/silent)
      for ix=0,31 do for iy=0,31 do begin
        ii=where(pixpos.ref_det1x ne -1 and pixpos.rawx eq ix and $
              pixpos.rawy eq iy and pixpos.grade le 26)

        thispdf=fltarr(360,360)

        if ii[0] ne -1 then for i=0,n_elements(ii)-1 do $
              if finite(total(pixpos[ii[i]].pdf)) then $
              thispdf[pixpos[ii[i]].ref_det1x:pixpos[ii[i]].ref_det1x+6,$
                    pixpos[ii[i]].ref_det1y:pixpos[ii[i]].ref_det1y+6]+= $
                    pixpos[ii[i]].pdf*gw[pixpos[ii[i]].grade]
        ii=where(thispdf gt allpdf)
        if ii[0] ne -1 then begin
            allpdf[ii]=thispdf[ii]
            pixmapx[ii]=ix
            pixmapy[ii]=iy
            detnum[ii]=idet
        endif

; Older code chunk, ignores grade
;    x=pixpos[ii].ref_det1x-1
;    y=pixpos[ii].ref_det1y-1
;    pdf=pixpos[ii].pdf
;    grade=pixpos[ii].grade
;    rawx=pixpos[ii].rawx
;    rawy=pixpos[ii].rawy
;    for i=0,n_elements(x)-1 do if finite(total(pdf[*,*,i])) then begin
;        jj=where(pdf[*,*,i] gt 0.0)
;;        jj2d=array_indices(intarr(7,7),jj)
;        thispdf=intarr(7,7)
;        thispdf[*,*]=-1
;        thisdet=thispdf
;        thispdf[jj]=rawx[i]+rawy[i]*32
;        thisdet[jj]=idet
;        pixmap[x[i]:x[i]+6,y[i]:y[i]+6]=thispdf
;        detnum[x[i]:x[i]+6,y[i]:y[i]+6]=thisdet
;    endif

      endfor
    endfor
    fits_write,auxildir+'pixmapx'+ab+'.fits',pixmapx,header
    fits_write,auxildir+'pixmapy'+ab+'.fits',pixmapy,header
    fits_write,auxildir+'detnum'+ab+'.fits',detnum,header
endif else begin
    fits_read,auxildir+'pixmapx'+ab+'.fits',pixmapx,header
    fits_read,auxildir+'pixmapy'+ab+'.fits',pixmapy,header
    fits_read,auxildir+'detnum'+ab+'.fits',detnum,header
endelse
end
