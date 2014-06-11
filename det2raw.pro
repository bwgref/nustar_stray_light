; converts from detector 64x64 array to RAWX/Y
pro det2raw,fp,fpd
npix = 64

fpd = lonarr(npix/2.,npix/2.,4) 
fpd[*,*,1]=fp[0:31,32:63]
fpd[*,*,0]=fp[32:63,32:63]
fpd[*,*,2]=fp[0:31,0:31]
fpd[*,*,3]=fp[32:63,0:31]

fpd[*,*,0] = rotate(reverse(fpd[*,*,0],2),1)  ; reverse y axis and then rotate by 90deg 
fpd[*,*,1] = reverse(fpd[*,*,1],1)            ; reverse x axis
fpd[*,*,2] = rotate(reverse(fpd[*,*,2],1),1)  ; reverse x axis and then rotate by 90deg 
fpd[*,*,3] = reverse(fpd[*,*,3],2)            ; reverse y axis

end

