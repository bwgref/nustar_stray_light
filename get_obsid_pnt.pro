;;
;; Finds obsid for a given position
;; Ex.: 
;;    get_obsid_pnt, 'SgrABridgePointing.txt', 266.460417 , -28.824444  , PA=196.0
;;
;; Pointings file should be in the following format: 
;; "obsid [string] ra [float] dec [float]"
;;

;;
;; ATTENTION! There is confusion with PA definition
;; PAsky = -180 -PAinst
;; 

pro get_obsid_pnt, file, ra, dec, PAreq=PAreq

  read_nustar_pointings, file,pnt_ra,pnt_dec,pnt_lon=pnt_lon,pnt_lat=pnt_lat, obsid=obsid
  
  PAsky=-180-PAreq
  
  for i=0,n_elements(obsid)-1 do begin
     make_astr_nustar, pnt_ra(i), pnt_dec(i), PAsky, astr=astr
     ad2xy, ra, dec, astr, x, y
     x=fix(x)
     y=fix(y)
     if(x ge 0 and x le 63 and y ge 0 and y le 63) then begin
        ret='yes'
     endif else begin 
        ret='no'
     endelse
     print,obsid(i),' ',ret
  endfor

end
