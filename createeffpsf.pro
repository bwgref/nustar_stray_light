;------------------------------------ 
; NAME: 
;	CreateEffPSF.pro 
;
; PURPOSE: 
;	This program creates and outputs a weighted PSF based on the off-axis histogram (as
;	outputted by nuproducts) at a given pixel location. It calls EffPSF.pro
;
; INPUTS:
;	x0 = x-location of the center of the region used to create the histogram [pixel]
;	y0 = y-location of the center of the region used to create the histogram [pixel]
;	PA = PA angle in degrees (from the PA_PNT header keyword in the NuSTAR fits file)
;	xi = x-location of the pixel for which the effective PSF is being created
;	yi = y-location of the pixel for which the effective PSF is being created
;	off_histo_name = name and path(if needed) of the off-axis histogram
;	
;	
; OUTPUTS:
;	PSFeff = 325x325 PSF image.  It has the same degPerPixel that the PSF uses.
;
; EXAMPLE:
;	(500, 499) = the central pixel for region used to create the off-axis histogram
;	(300, 321) = the pixel for which you want to create the effective PSF
;	PA = 0 deg
;	IDL> CreateEffPSF, 500, 499, 0, 300, 321, 'offaxishistB.fits', PSFeff = myPSF
;
; CAUTION: 
;	The path for the PSF has been hard-coded in below.  Please change it to your own path.
;
; MODIFICATIONS:
;	May 2, 2013:  initial release
;
; AUTHOR:
;	Melania Nynka mcd2005@columbia.edu
;			
;------------------------------------------------------------------------------


PRO CreateEffPSF, x0, y0, PA, xi, yi, off_histo_name, PSFeff=PSFeff, hdr_psf=hdr_psf, $
                  astr_psf=astr_psf, hdr_sky=hdr_sky, arcminPerPixPSF=arcminPerPixPSF, psf_koeff=psf_koeff

if(getenv('CALDB') eq '') then begin
   print,'*** ERROR ***'
   print,'Hey, where is NuSTAR $CALDB variable?'
   stop
endif else begin
   print,'Using $CALDB=',getenv('CALDB')
endelse


PSFpath = getenv('CALDB')+'/data/nustar/fpm/bcf/psf/nuA2dpsf20100101v003.fits'

OffAxisAngles=indgen(18)/2.0		
dr = !PI/180.
rd = 180./!PI

xyad,hdr_sky,xi,yi,xi_ra,xi_dec



  if(n_elements(psf_koeff) eq 0) then psf_koeff=1.0      

;;;;;;;;;;;; load  the 18 PSFs

;; creates a psf array 
increment = (OffAxisAngles[1]-OffAxisAngles[0])

;; loads the images and headers into arrays
for i=0, (size(OffAxisAngles))[1]-1 do begin	

   if(~file_test(PSFpath)) then begin
      print, 'CALDB file does not exist: ',PSFpath
      stop
   endif

	tempPSF= mrdfits(PSFpath, i+1, hdr, /silent)
	if i eq 0 then begin 
		hdrArray = hdr
		psfSizeX = fxpar(hdr, 'NAXIS1')
		psfSizeY = fxpar(hdr, 'NAXIS2')
		arcminPerPixPSF = 60.0*fxpar(hdr, 'CDELT1')  ;symmetric x,y
		psfArray = make_array(psfSizeX, psfSizeY, (size(OffAxisAngles))[1], value=0.0)
	endif
	
	tempPSF = tempPSF/total(tempPSF)  ; normalize the psf to total=1
	psfArray[*,*,i] = tempPSF	
	
endfor


;;;;;;;;;;;;; load and filter off-axis histos and trim zeros

if (~file_test(off_histo_name)) then begin
   print, 'Not found: ',off_histo_name
   stop
endif

histos = mrdfits(off_histo_name, 'OFFAXIS_HISTO', histhear,/SILENT)
theta = (histos.off_axis)[where(histos.duration ne 0.0)]	; arcmins (radius)
phi = (histos.phi)[where(histos.duration ne 0.0)]			; degrees (rotation angle)
time = (histos.duration)[where(histos.duration ne 0.0)]		; seconds (exposure time)



;; find the x,y shifts and average onaxis4 positions (convert theta->pixels, phi->radians)
totRot = phi+(360-PA)
xOnaxisShift = (theta/arcminPerPixPSF)*cos(totRot*3.14159/180.0)	; in pixels
yOnaxisShift = (theta/arcminPerPixPSF)*sin(totRot*3.14159/180.0)	; in pixels
xOnaxisPos = x0-xOnaxisShift                                            ; pixels							
yOnaxisPos = y0-yOnaxisShift                                            ; pixels



EffPSF, psfArray, OffAxisAngles, arcminPerPixPSF, xOnaxisPos, yOnaxisPos, time, xi, yi, ePSF=PSFeff
sz=size(PSFeff)
sx=sz[1]
sy=sz[2]
alpha=0
delt=arcminPerPixPSF/60.*psf_koeff
crpix=[sx/2.,sy/2.]
crval=[xi_ra,xi_dec]
rot_mat = [ [ cos(alpha*dr), -sin(alpha*dr)], $ ;Rotation matrix
            [ sin(alpha*dr),  cos(alpha*dr)] ] 
make_astr,astr_psf, CD=rot_mat, DELT = [-delt, delt], CRPIX = crpix, CRVAL = crval,$
          RADECSYS = 'FK5', EQUINOX = 2000.0
mkhdr, hdr_psf, 4, [sx,sy]
putast, hdr_psf, astr_psf, CD_TYPE=2
SXADDPAR, hdr_psf, 'PNT_PA', PA, ' PA angle, deg.'
SXADDPAR, hdr_psf, 'ALPHA', alpha, ' Alpha angle, deg.'
SXADDPAR, hdr_psf, 'RA_OBJ', xi_ra, ' RA'
SXADDPAR, hdr_psf, 'DEC_OBJ', xi_dec, ' Dec'

print, "finished creating effective PSF"


end
