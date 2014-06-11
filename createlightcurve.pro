;+
;NAME:
;
; createlightcurve.pro
; 
; 
;PURPOSE:
;
; Create a regulatrly binned lightcurve for either one or both foal planes 
; of NuSTAR given the event files. There are options to correct for livetime 
; and vigneting factor. It is also possible to produce a lightcurve where 
; each bin is a GTI (usually an orbit).
; 
; 
;CALL:
;
; createlightcurve( evtfn[0], binsize=150, Elo=Emin, Ehi=Emax, $
;                      livetimefn=livetimefn[0], offaxisfn=offaxisfn[0], photindex=2.3,  $
;                      _eventfn=evtfn[1], $
;                      _livetimefn=livetimefn[1], _offaxisfn=offaxisfn[1], _photindex=2.3 )
;
;
;INPUT:
;
; eventfn     : event filename 
; 
;OPTIONAL INPUT:
;
; outpath     : output directory (default is ./)                   
; binsize     : [s] size of the bins (default is 100s)
; Elo         : [keV] Low boundary of energy range (default is 3 keV)
; Ehi         : [keV] High boundary of energy range (default is 79 keV)
; livetimefn  : File name and path of the fits file containing the livetime info 
; offaxisfn   : File name and path of the fits file containing the offaxis info
; photindex   : Photon index of the source for vigneting correction
; GTIaverage  : Switch to enable lightcurve with irregular bins matching the GTIs (one bin per GTI)
; _eventfn    : Same for second focal plane
; _livetimefn :
; _offaxisfn
; _photindex
;
;
;OUTPUT:
; 
;
;  
;  
;EXTERNAL CALLS:
;
; FXBOPEN  (part of GSFC's astron library)
; FXBREAD  (part of GSFC's astron library)
; FXBCLOSE (part of GSFC's astron library)
; HEADFITS (part of GSFC's astron library)
; FXPAR  (part of GSFC's astron library)
; upper_limit_gehrels  (Nicolas Barriere's )
; lower_limit_gehrels  (Nicolas Barriere's )
; GetBinnedLivetime  (Nicolas Barriere's )
; GetBinnedVignet  (Nicolas Barriere's )
;
;
;HISTORY:
;Roman Krivonos -- vigntting was commented out with ";=;",
;                  krivonos@ssl.berkeley.edu Feb 2014 
;Nicolas Barrière, University of California at Berkeley, barriere@ssl.berkeley.edu
;October 2012
;-



FUNCTION Createlightcurve, eventfn, outpath=outpath, $
                      binsize=binsize, Elo=Elo, Ehi=Ehi, $
                      livetimefn=livetimefn, offaxisfn=offaxisfn, photindex=photindex,  $
                      _eventfn=_eventfn, $
                      _livetimefn=_livetimefn, _offaxisfn=_offaxisfn, _photindex=_photindex, $
                      GTIaverage=GTIaverage, $
                      detector=detector

IF ~KEYWORD_SET(binsize) THEN binsize=100.
IF ~KEYWORD_SET(Elo) THEN Elo=3.
IF ~KEYWORD_SET(Ehi) THEN Ehi=79.
IF ~KEYWORD_SET(phindex) THEN phindex = 2.
IF ~KEYWORD_SET(phindex) THEN phindex = 2.
IF ~KEYWORD_SET(outpath) THEN outpath = './'
IF ~KEYWORD_SET(GTIaverage) THEN GTIaverage = 0


;===============================================================
;   Read input file, load event times and start and stop times
;===============================================================
FXBOPEN, lun, eventfn, 1
FXBREAD, lun, tt, 'TIME'  ;[s]  Event Time (seconds since Jan 2010 00:00:00 UTC)
FXBREAD, lun, pi, 'PI'      ;[channel] Event Pulse Invariant
FXBREAD, lun, detID, 'DET_ID'   
FXBCLOSE, lun


FXBOPEN, lun, eventfn, 2   ; GTI extension
FXBREAD, lun, GTIstart, 'START' 
FXBREAD, lun, GTIstop, 'STOP'      
FXBCLOSE, lun
nGTI = N_ELEMENTS(GTIstart)


header = HEADFITS(eventfn)
module = STRMID(FXPAR(header, 'INSTRUME'), 3,1)
IF ~KEYWORD_SET(tstart) THEN tstart = FXPAR(header, 'TSTART')
IF ~KEYWORD_SET(tstop) THEN tstop = FXPAR(header, 'TSTOP')

IF KEYWORD_SET(_eventfn) THEN BEGIN
    FXBOPEN, lun, _eventfn, 1
    FXBREAD, lun, tt2, 'TIME'  ;[s]  Event Time (seconds since Jan 2010 00:00:00 UTC)
    FXBREAD, lun, pi2, 'PI'      ;[channel] Event Pulse Invariant
    FXBCLOSE, lun
    
    header = HEADFITS(_eventfn)
    _module = STRMID(FXPAR(header, 'INSTRUME'), 3,1)
ENDIF



;===============================================================
;   Filter energy
;===============================================================
en = 0.04*pi+1.6                         ;[keV] events energy

igood = WHERE((en GE Elo) AND (en LE Ehi), ngood)
tt = tt[igood]
detID = detID[igood]

IF n_elements(detector) NE 0 THEN BEGIN
  igood = where(detID EQ detector)
  tt = tt[igood]
ENDIF

IF KEYWORD_SET(_eventfn) THEN BEGIN
    en2 = 0.04 * pi2 + 1.6
    igood2 = WHERE((en2 GE Elo) AND (en2 LE Ehi), ngood)
    tt2 = tt2[igood2]
ENDIF


;===============================================================
;   Build lightcurve, get upper and lower limits
;===============================================================
IF ~GTIaverage THEN BEGIN
    rate = HISTOGRAM(tt, BINSIZE=binsize, MIN=tstart, MAX=tstop, LOCATIONS=binstart) 
    binsize = REPLICATE(binsize, N_ELEMENTS(rate))
ENDIF ELSE BEGIN
    rate = DBLARR(nGTI)
    FOR i=0, nGTI-1 DO BEGIN

        igood = WHERE((tt GT GTIstart[i]) AND (tt LE GTIstop[i]), ngood)
        rate[i] = ngood

    ENDFOR
    binstart = GTIstart
    binsize = GTIstop - GTIstart
ENDELSE


IF KEYWORD_SET(_eventfn) THEN BEGIN
    rate2 = HISTOGRAM(tt2, BINSIZE=binsize[0], MIN=tstart, MAX=tstop)  
    
    upperlim = upper_limit_gehrels(0.8413, rate+rate2)
    lowerlim = lower_limit_gehrels(0.8413, rate+rate2)

ENDIF ELSE BEGIN

    upperlim = upper_limit_gehrels(0.8413, rate)
    lowerlim = lower_limit_gehrels(0.8413, rate)
    
ENDELSE

;==============================================================
;   Remove LC bins that are out of GTIs
;==============================================================
FOR i=0, nGTI-2 DO BEGIN

    iout = where(binstart GE GTIstop[i] and binstart+binsize LE GTIstart[i+1], nout, COMPLEMENT=iin)
    binstart = binstart[iin]
    binsize = binsize[iin]
    rate = rate[iin]
    IF KEYWORD_SET(_eventfn) THEN rate2 = rate2[iin]
    upperlim = upperlim[iin]
    lowerlim = lowerlim[iin]
    
ENDFOR


;==============================================================
;   Correct for bins that are partially out of GTIs
;===============================================================
FOR i=0, nGTI-1 DO BEGIN

    igtistart = WHERE(GTIstart[i] GT binstart AND GTIstart[i] LT (binstart+binsize), ngtistart)
    IF ngtistart NE 0 THEN BEGIN
        binsize[igtistart] = binsize - (GTIstart[i] - binstart[igtistart])
        binstart[igtistart] = GTIstart[i]
    ENDIF
    
    
    igtistop = WHERE(GTIstop[i] GT binstart AND GTIstop[i] LT (binstart+binsize), ngtistop)
    IF ngtistop NE 0 THEN BEGIN
        binsize[igtistop] = GTIstop[i] - binstart[igtistop]
    ENDIF
ENDFOR


itooshort = WHERE(binsize LT 0.2*MEAN(binsize), ntooshort, COMPLEMENT=iok)
binstart = binstart[iok]
binsize = binsize[iok]
rate = rate[iok]
IF KEYWORD_SET(_eventfn) THEN rate2 = rate2[iok]
upperlim = upperlim[iok]
lowerlim = lowerlim[iok]


;===============================================================
;   Livetime and vigneting factor averaged for each lightcurve bin
;===============================================================
IF KEYWORD_SET(livetimefn) THEN BEGIN
    livetime = GetBinnedLivetime( livetimefn, binstart, binsize )
    rate /= livetime
    livestring = '_livecorr'
    
    IF KEYWORD_SET(_livetimefn) THEN BEGIN
        livetime2 = GetBinnedLivetime( _livetimefn, binstart, binsize )
        rate2 /= livetime2
        upperlim /= (livetime + livetime2)/2.
        lowerlim /= (livetime + livetime2)/2.
    ENDIF ELSE BEGIN
        upperlim /= livetime
        lowerlim /= livetime
    ENDELSE
     
ENDIF ELSE livestring = ''
    

IF KEYWORD_SET(offaxisfn) THEN BEGIN
;=; vignet = GetBinnedVignet(offaxisfn, binstart, binsize, Elo, Ehi, photindex=photindex, FPM=module)
   print,'Ask Nicolas for GetBinnedVignet function'
   stop
   rate /= vignet
   vignetstring = '_vigncorr'
  
   IF KEYWORD_SET(_offaxisfn) THEN BEGIN
;=; vignet2 = GetBinnedVignet(_offaxisfn, binstart, binsize, Elo, Ehi, photindex=photindex, FPM=_module)
   print,'Ask Nicolas for GetBinnedVignet function'
   stop
   rate2 /= vignet2
   upperlim /= (vignet + vignet2) / 2.
   lowerlim /= (vignet + vignet2) / 2.
  ENDIF ELSE BEGIN
     upperlim /= vignet
     lowerlim /= vignet
  ENDELSE
  
ENDIF ELSE BEGIN
  vignetstring = ''
  
ENDELSE
  

;===============================================================
;   ; lightcurve in [cts/s]  !!! corrected counts, not anymore raw counts if livetime and vigneting corrections applied
;===============================================================
upperlim /= binsize
lowerlim /= binsize

IF KEYWORD_SET(_eventfn) THEN rate = (rate + rate2) / binsize  ELSE rate /= binsize 
  
 
;===============================================================
; Look for NANs and set them to 0.
;===============================================================
inan = WHERE(FINITE(rate) EQ 0, nnan)
IF nnan NE 0 THEN rate[inan] = 0.

inan = WHERE(FINITE(upperlim) EQ 0, nnan) 
IF nnan NE 0 THEN upperlim[inan] = 0.

inan = WHERE(FINITE(lowerlim) EQ 0, nnan) 
IF nnan NE 0 THEN lowerlim[inan] = 0.


;===============================================================
;   Store result in a sav file
;===============================================================
lc = {rate      :rate,      $
      upperlim  :upperlim,  $
      lowerlim  :lowerlim,  $
      tstart    :tstart,    $
      tstop     :tstop,     $  
      binstart  :binstart,  $
      binsize   :binsize,   $
      Elo       :Elo,       $
      Ehi       :Ehi        }       



RETURN, lc

END
