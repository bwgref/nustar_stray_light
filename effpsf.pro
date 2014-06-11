;------------------------------------ 
; NAME: 
;	EffPSF.pro 
;
; PURPOSE: 
;	This program creates and outputs a weighted PSF based on the off-axis histogram (as
;	outputted by nuproducts) at a given pixel location.
;
; INPUTS:
;	psfArray = an array of PSF images.  size = [PSFx, PSFy, #angles]
;				See CreateEffPSF or Convolution2 for an example 
;	OffAxisAngles = an array of PSF off-axis angles.  
;	aminPerPix = arcminute p0er pixel. CDELT1*60.0 header keyword from the PSF
;	xoa_array = array of x-onaxis positions.  See CreateEffPSF or Convolution2 for an example 
;	yoa_array = array of x-onaxis positions.  See CreateEffPSF or Convolution2 for an example
;	time_array = array of exposure time for the x,y onaxis positions.  
;			See CreateEffPSF or Convolution2 for an example
;	xi = x pixel position for which you want to create the effective PSF
;	yi = y pixel position for which you want to create the effective PSF
;	
;	
; OUTPUTS:
;	ePSF = 325x325 PSF image.  It has the same degPerPixel that the PSF uses.
;
; EXAMPLE:
;	(500, 499) = the central pixel for region used to create the off-axis histogram
;	(300, 321) = the pixel for which you want to create the effective PSF
;	PA = 0 deg
;	IDL> CreateEffPSF, 500, 499, 0, 300, 321, 'offaxishistB.fits', PSFeff = myPSF
;
;
; MODIFICATIONS:
;	May 2, 2013:  initial release
;
; AUTHOR:
;	Melania Nynka mcd2005@columbia.edu
;			
;------------------------------------------------------------------------------

Pro EffPSF, psfArray, OffAxisAngles, aminPerPix, xoa_array, yoa_array, time_array, xi, yi, ePSF=ePSF

;; create a psf-sized image with all pixels set to 0.
ePSF = make_array( (size(psfArray[*,*,1]))[1], (size(psfArray[*,*,1]))[1], value=0.0)	
	
	
;; check of on-axis x,y,time arrays are the same size
if ( ((size(xoa_array))[1] ne (size(yoa_array))[1]) || $
	((size(xoa_array))[1] ne (size(time_array))[1]) || $
	((size(yoa_array))[1] ne (size(time_array))[1]) ) then begin
	print, "ERROR: The on-axis arrays are not the same size."
	print, "x-onaxis array: ", size(xoa_array)
	print, "y-onaxis array: ", size(yoa_array)
	print, "time array: ", size(time_array)
endif


for i=0, (size(xoa_array))[1]-1 do begin
	
	xoa = xoa_array[i]	; pixels
	yoa = yoa_array[i]	; pixels
	time = time_array[i]	

	delx = xi-xoa
	dely = yi-xoa
	increment = (OffAxisAngles[1]-OffAxisAngles[0])
		
	theta_ij = sqrt((delx)^2+(dely)^2)			; radius in pixels
	thetaAmin_ij = theta_ij*aminPerPix			; radius in arcmin
	phi_ij = (atan(dely/delx))*(180.0/3.14159)	; CW rotation in deg
	
	;; choosing the PSF image
	floorRadius = max(where(OffAxisAngles le thetaAmin_ij))  ; defines index for off-axis-angles
	if (floorRadius lt (size(OffAxisAngles))[1]-1 && (thetaAmin_ij-OffAxisAngles[floorRadius]) gt increment/2.) then begin
		floorRadius++
	endif
	
	PSFimg = psfArray[*,*,floorRadius]
	
	if phi_ij ne 0.0 then PSFimg = ROT(PSFimg, phi_ij)
	
	ePSF = ePSF+(PSFimg*time)
	
	;;print, thetaAmin_ij, OffAxisAngles[floorRadius], phi_ij, time, xoa, yoa

endfor


end
