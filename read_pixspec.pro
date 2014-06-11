pro read_pixspec, path=path, e1=e1, e2=e2, fpA=fpA, fpB=fpB

  if(n_elements(e1) eq 0) then e1=3.
  if(n_elements(e2) eq 0) then e2=80.

  cha1=LONG((e1-1.6)/0.04)
  cha2=LONG((e2-1.6)/0.04)
  
  

  det1w=360
  npix = 64

  fpdA = lonarr(npix/2,npix/2,4) 
  fpdB = lonarr(npix/2,npix/2,4) 

  module=['A','B']
  for imod=0,1 do begin
     for id=0,3 do begin
        for x=1,30 do begin
           for y=1,30 do begin
              file=path+'/PixSpec_'+module[imod]+'_'+String(id,Format='(i02)')+'_'+String(x,Format='(i02)')+'_'+String(y,Format='(i02)')+'.fits'
              if(~file_test(file)) then begin
                 print,'file? ',file
                 stop
              endif
              spec = mrdfits(file, 1, /FSCALE, header,/SILENT)
              if(imod eq 0) then fpdA[x,y,id]=TOTAL((spec.normcounts)[(cha1-1):(cha2-1)])
              if(imod eq 1) then fpdB[x,y,id]=TOTAL((spec.normcounts)[(cha1-1):(cha2-1)])
           endfor
        endfor
     endfor
  endfor

  raw2det,fpdA,fpA
  raw2det,fpdB,fpB

end
