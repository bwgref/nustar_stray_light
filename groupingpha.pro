;+
;NAME:
;
; groupingPHA.pro
; 
; 
;PURPOSE:
;
; Groups spectral bins by statistical significance.
; Produces a ASCII file that is then used by grppha to create the actual grouped spectrum fits file
; 
; 
;CALL:
;
; groupingPHA, 'source.pha', 'background.pha', 'groups.dat', nsigma = 2., chanmax=960
;
;
;INPUT:
;
; inputSRC  : PHA filename for the source
; inputBKG  : PHA filename for the background
; outfile   : Name of the ASCII file that contains the grouping (to be input in grppha) 
; 
; 
;OPTIONAL INPUT:
;
; path          : path to the PHA file (default is ./)
; nsigma        : Significance level in each group of bins (default is 3)           
; chanmin       : lowest chanel to take into account (inclusive) 
; chanmax       : Highest channel to take into account (inclusive) 
; lastbinsigmin : Minimum significance level in the last group. If this 
;               : level is not reached, then this bin is merged with the previous one.
;
;
;OUTPUT:
; 
; Writtes the groups in a n ASCII file (filename provided in call). Also print a bunch of info in the 
; terminal, including one line that should be copied and pasted in the terminal to execute grppha.
;  
;
;HISTORY:
;
; Loosely based on the code of the same name written by Melania Nynka and Kerstin Perez (Columbia University)
;
; Nicolas Barrere, UC Berkeley, Space Sciences laboratory, barriere@ssl.berkeley.edu
; November 2012
;
;-


PRO groupingPHA, inputSRC, inputBKG, outfile, path=path, nsigma=nsigma, $
                 chanmin=chanmin, chanmax=chanmax, lastbinsigmin=lastbinsigmin, $
                 skipback=skipback, backnorm=backnorm, silent=silent, exec=exec, exposure_key=exposure_key

  if(n_elements(silent) eq 0) then silent=0
  if(n_elements(exec) eq 0) then exec=0 else exec=1
  if(n_elements(exposure_key) eq 0) then exposure_key='LIVETIME'

;===============================
;  Default values for optional parameters
;===============================
IF ~KEYWORD_SET(skipback) THEN skipback = 0 ELSE skipback = 1
IF ~KEYWORD_SET(path) THEN path = './'
IF ~KEYWORD_SET(nsigma) THEN nsigma = 3.
IF ~KEYWORD_SET(chanmin) THEN chanmin = 35      ; 3 keV for NuSTAR
IF ~KEYWORD_SET(chanmax) THEN chanmax = 1935    ; 79 keV for NuSTAR


;===============================
;  Read fits files
;===============================
cd, path, CURRENT=currentdir

specSRC = mrdfits(inputSRC, 1, /FSCALE, srcheader)
specBKG = mrdfits(inputBKG, 1, /FSCALE, bkgheader)
print,Total(specBKG.counts)


backscalBKG = fxpar(bkgheader, 'BACKSCAL')
backscalSRC = fxpar(srcheader, 'BACKSCAL')

areascalBKG = fxpar(bkgheader, 'AREASCAL')
areascalSRC = fxpar(srcheader, 'AREASCAL')

livetimeBKG = fxpar(bkgheader, exposure_key)
livetimeSRC = fxpar(srcheader, exposure_key)


;===============================
;  Net counts
;===============================
print,backscalSRC,livetimeSRC,areascalSRC,' / ', BackscalBKG , livetimeBKG , areascalBKG
backnorm = ( backscalSRC * livetimeSRC * areascalSRC / (backscalBKG * livetimeBKG * areascalBKG) )


if(skipback eq 1) then begin
   backnorm=0.0
   if(silent ne 1) then begin
      print,'***'
      print,'*** WARNING: ignoring background (backnorm=0)'
      print,'***'
   endif
endif

subtracted = specSRC.counts - specBKG.counts * backnorm ;> 0

;Total net counts
totalnet = TOTAL(subtracted[chanmin:chanmax])
totalerror = SQRT( TOTAL((specSRC.counts)[chanmin:chanmax]) + TOTAL((specBKG.counts)[chanmin:chanmax]) * backnorm^2 )
totallevel = totalnet / totalerror


;===============================
;  Initializations
;===============================
ChannelSize=(size(subtracted))[1]
tempsum=0.0
minchan = [0]
maxchan = [0]
minchanholder = 0
counter = 1


groupstart = 0
groupstop = 0
grouplevel = 0D
groupcounts = 0D
grouperror = 0D


;===============================
;  Group the bins
;===============================
first = 1
ngroups = 0

REPEAT BEGIN

  IF first THEN BEGIN
    istart = chanmin 
    first = 0
  ENDIF ELSE BEGIN
    istart = istop+1
  ENDELSE
  
  istop = istart
  
  REPEAT BEGIN
    signalcounts = TOTAL(subtracted[istart:istop])
    error = SQRT( TOTAL((specSRC.counts)[istart:istop]) + TOTAL((specBKG.counts)[istart:istop]) * backnorm^2 )
    level = signalcounts / error
    
    IF istop EQ chanmax THEN BREAK ELSE istop++
  ENDREP UNTIL level GE nsigma
  
  groupstop = [groupstop, istop]
  groupstart = [groupstart, istart]
  grouplevel = [grouplevel, level]
  groupcounts = [groupcounts, signalcounts]
  grouperror = [grouperror, error]
  
  ngroups++
  
ENDREP UNTIL istop EQ chanmax
  
  
;===============================
;  Remove first bin of each array (was there only to initialize them)
;===============================
groupstart = groupstart[1:ngroups]
groupstop = groupstop[1:ngroups]
grouplevel = grouplevel[1:ngroups]
groupcounts = groupcounts[1:ngroups]
grouperror = grouperror[1:ngroups]


;===============================
;  Check last bin
;===============================
IF KEYWORD_SET(lastbinsigmin) THEN BEGIN
  IF grouplevel[ngroups-1] LT lastbinsigmin THEN BEGIN
    ; Append the two last groups in the second to last, and delete the last group
  
    groupstart = groupstart[0:ngroups-2]
    groupstop[ngroups-2] = groupstop[ngroups-1]
    groupstop = groupstop[0:ngroups-2]
    
    groupcounts[ngroups-2] = TOTAL(groupcounts[ngroups-2:ngroups-1])
    groupcounts = groupcounts[0:ngroups-2]
    
    grouperror[ngroups-2] = SQRT( TOTAL((specSRC.counts)[groupstart[ngroups-2]:groupstop[ngroups-2]]) + $
            TOTAL((specBKG.counts)[groupstart[ngroups-2]:groupstop[ngroups-2]]) * backnorm^2 )
    grouperror=grouperror[0:ngroups-2]
    
    grouplevel[ngroups-2] = groupcounts[ngroups-2] / grouperror[ngroups-2]
    grouplevel = grouplevel[0:ngroups-2]
    
    ngroups = ngroups-1
    
  ENDIF

ENDIF


;===============================
;  Print ASCII file 
;===============================
OPENW, lun, outfile, /get_lun  

FOR j=0, ngroups-1 DO BEGIN
  difference = groupstop[j]-groupstart[j]+1
  PRINTF, lun, groupstart[j], groupstop[j], difference, FORMAT='(I4, 4X, I4, 4X, I4)'
ENDFOR

FREE_LUN, lun

if(silent ne 1) then begin
;===============================
;  Print results in terminal
;===============================
   PRINT, '**************************************************'
   PRINT, 'start, stop, n chan, sigma, signal, error' 
   FOR j=0, ngroups-1 DO PRINT, groupstart[j], groupstop[j], groupstop[j]-groupstart[j]+1, $
                                grouplevel[j], groupcounts[j], grouperror[j], $
                                format='(I4, 2X, I4, 4X, I4, 3X, F5.2, 1X, F6.1, 2X, F5.1)'

   PRINT, '**************************************************'
endif

PRINT, totalnet, totalerror, totallevel, FORMAT='("Total net counts: ", F10.1, " +/- ", F10.1, " (", F10.1, " sigma)" )'

if(silent ne 1) then begin
   PRINT, ''
   PRINT, "*** groupingPHA: wrote grouping file " + outfile
endif

grpfn = STRSPLIT(inputSRC , '/.', /EXTRACT) 
grpfn = grpfn[n_elements(grpfn)-2] + '_grp.fits'

PRINT, ''
PRINT, 'You may want to execute this line:'
cmd = 'grppha ' +  inputSRC + ' ' + grpfn + ' clobber=yes comm="bad 0-' + STRTRIM(STRING(chanmin, FORMAT= '(I10)'), 2) + $
      ' & group ' +  outfile + ' & bad ' + STRTRIM(STRING(chanmax, FORMAT= '(I10)'), 2) + '-4095 & good ' + $
      STRTRIM(STRING(chanmin, FORMAT= '(I10)'), 2) + '-' + STRTRIM(STRING(chanmax, FORMAT= '(I10)'), 2)+' & exit" '
PRINT, cmd
PRINT, ''

if(exec eq 1) then spawn,cmd

CD, currentdir

END







