; converts from RAWX/Y to 64x64 detector array 
pro raw2det,fpd,fp
npix = 64
type=Size(fpd,/TYPE)
fp = make_array(npix,npix,type=type)	; One complete focal plane 
fptemp = make_array(npix/2.,npix/2.,4,type=type) 

fptemp[*,*,0] = reverse(rotate(fpd[*,*,0],3),2)  ; rotate by 270deg  and then reverse y axis
fptemp[*,*,1] = reverse(fpd[*,*,1],1)            ; reverse x axis
fptemp[*,*,2] = reverse(rotate(fpd[*,*,2],3),1)  ; rotate by 270 deg and then reverse x axis
fptemp[*,*,3] = reverse(fpd[*,*,3],2)            ; reverse y axis

fp[0:31,32:63]  = fptemp[*,*,1]          ;top left corner
fp[32:63,32:63] = fptemp[*,*,0]          ;top right corner
fp[0:31,0:31]   = fptemp[*,*,2]          ;bottom left
fp[32:63,0:31]  = fptemp[*,*,3]         ; bottom right

end

