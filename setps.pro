PRO SETPS, portrait = portrait
; Sets the plotting environment for PS plots:
  COMMON PS_BLOCK, outps, outpdf

  ; Get the current working directory
  cd, current = dir

  IF  ~keyword_set(outps) THEN BEGIN
     outps=dir+'/idlout.ps'
     outpdf = dir+'/idlout.pdf'
  ENDIF

  IF ~keyword_set(portrait) THEN BEGIN
     xsize = 9
     ysize = 6
     xoffset = 1.25
     yoffset = 10
     landscape = 1
     !p.charsize = 1
  ENDIF ELSE BEGIN
     xsize = 8
     ysize = 10
     xoffset = 0.25
     yoffset = 0.5
     landscape = 0
     !p.charsize = 1.5
  ENDELSE


  set_ps_plot=1
  set_plot,'ps'
  !p.thick = 3
  !p.charthick = 3
  !p.multi = 0

  device,filename=outps,/inches,/color,landscape = landscape,$
         xsize=xsize,ysize=ysize,xoffset=xoffset,yoffset=yoffset
  
end
