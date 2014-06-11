pro print_gc_grid_eq
  deg = 3.14159265359/180
  arcmin = deg/60.0

  HalfField = 6.7*arcmin

  ra_max=268.42
  ra_min=263.93
  dec_min=-31.2
  dec_max=-27.47

  LatMin = -0.4*deg + HalfField
  LatMax = +0.4*deg - HalfField

  LongMin = -1*deg + HalfField
  LongMax = +1*deg - HalfField
  
  LatPointings = floor((LatMax-LatMin)/HalfField)+1
  LongPointings = floor((LongMax-LongMin)/HalfField)+1


  if(LatPointings MOD 2 eq 1) then LatPointings++
  if(LongPointings MOD 2 eq 1) then LongPointings++
  
  LatStart = -double(LatPointings-1)/2.0*HalfField
  LongStart = -double(LongPointings-1)/2.0*HalfField

  cnt=1
  openw, lun,  'print_gc_grid_eq.dat', /get_lun
  for ra=ra_min,ra_max,HalfField/deg do begin
     for dec=dec_min,dec_max,HalfField/deg do begin
        EULER, ra, dec, lon, lat, 1
        if(lon gt LongMin/deg and lon lt LongMax/deg and lat gt LatMin/deg and lat lt LatMax/deg) then begin
        ;EULER, lon, lat, ra, dec, 2
     
        printf,lun,'00000',cnt,'001',ra,dec,format='(a,i02,a,2f14.6)'
        cnt++
        endif
     endfor
  endfor
  close, lun
  free_lun, lun
  
end
