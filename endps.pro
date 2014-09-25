PRO ENDPS
; routine that closes the current postscript file and then spwans
; ps2pdf, which converts the PS file to a PDF.
  COMMON PS_BLOCK, outps, outpdf


  device,/close
  set_plot,'x'
  set_ps_plot=0
;  print,'/users/bwgref/bin/ps2pdf_idl.sh '+outps+' '+outpdf
  spawn,'./ps2pdf_idl.sh '+outps+' '+outpdf
  spawn, 'rm '+outps
  !p.multi=0
  
end
