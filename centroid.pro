;
; i,j - position of the maximum
;

; Gaussian sigma (assumed PSF)
smo=5./60./dabs(cdelt1)

;small step
dw=0.05
pmax=0
xwmax=0
ywmax=0
; some characteristic size
nn=2 

for xw=double(i-2),double(i+2),dw do begin
   for yw=double(j-2),double(j+2),dw do begin
      pcur=0
      do ii=i-nn,i+nn do begin
         do jj=j-nn,j+nn do begin
            ; check borders
            if(ii.ge.1.and.ii.le.nsx.and.jj.ge.1.and.jj.le.nsy) then begin
               xp=ii
               yp=jj
               weightx=dexp(-(xw-xp)**2/(2d0*smo**2))
               weighty=dexp(-(yw-yp)**2/(2d0*smo**2))
               pcur=pcur+map(ii,jj)*weightx*weighty
            endif
         endfor                     
      endfor                      
      
      if(pcur.gt.pmax) then begin
         pmax=pcur
         xwmax=xw
         ywmax=yw
      endif


   endfor                    
endfor                      
