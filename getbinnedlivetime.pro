;NAME:
;
; GetBinnedLivetime.pro
; 
; 
;PURPOSE:
;
; Get the livetime from the housekeeping file and rebin it to match the bins of the ligthcurve.
; 
; 
;CALL:
;
; result = GetEventLivetime(filename, binstart, binsize)
;
;
;INPUT:
;
; filename  : complete filename of the the FITS file containing the livetime
; binstart  : [s] start times of the lightcurve bins
; binsize   : [s] duration of the bins
; 
;
;OUTPUT:
; 
; This function returns an array containing as many bins as the binstart does.
; It returns the mean livetime over the each bin's duration.
; 
; !!!!  Note:
; Very dirty way to rebin in terms of positions of the boundaries, but accurate in terms of sum and pretty fast.
; It is OK as the bins are alsways larger than 100s
;  
;  
;EXTERNAL CALLS:
;
; fxbopen.pro   (part of the astro lib)
; fxbread.pro   (part of the astro lib)
; fxbclose.pro  (part of the astro lib)
; HSI_rebinner
;
;
;HISTORY:
;
; Nicolas Barriere, 10/04/2012
;
;-


FUNCTION GetBinnedLivetime, filename, binstart, binsize

FXBOPEN, lun, filename, 1
FXBREAD, lun, time, 'TIME'  ;[s]  FPM HK time column (secs since Jan 2010 00:00:00). The value in TIME marks the beginning of the 1-second period during w
FXBREAD, lun, livetime, 'LIVETIME'  ;  Shield corrected livetime for current second
FXBCLOSE, lun

livetimerebin = DBLARR(N_ELEMENTS(binstart))

FOR i=0, N_ELEMENTS(binsize)-1 DO BEGIN
  blah = MIN(ABS(binstart[i]-time), ilo)
  blah = MIN(ABS(binstart[i]+binsize[i]-time), ihi)
  livetimerebin[i] = TOTAL(livetime[ilo:ihi]) / FLOAT(ihi-ilo+1)
ENDFOR

RETURN, livetimerebin
END