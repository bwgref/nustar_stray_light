pro print_gc_grid
  deg = 3.14159265359/180
  arcmin = deg/60.0

  ;HalfField = 6.7*arcmin ; - original from Andreas
  ;HalfField = (6.0-2.26)*arcmin ; used for GC survey
  HalfField = (6.0-2.26)*arcmin*10

  LatMin = -0.4*deg + HalfField 
  LatMax = +0.4*deg - HalfField 
  LongMin = -1.17*deg + HalfField
  LongMax = +0.35*deg - HalfField
  LatPointings = floor((LatMax-LatMin)/HalfField)+1
  LongPointings = floor((LongMax-LongMin)/HalfField)+1

; test wide survey:
  LatMin = -5.*deg + HalfField 
  LatMax = +5.*deg - HalfField 
  LongMin = -10.*deg + HalfField
  LongMax = +10.*deg - HalfField
  LatPointings = floor((LatMax-LatMin)/HalfField)+1
  LongPointings = floor((LongMax-LongMin)/HalfField)+1


  if(LatPointings MOD 2 eq 1) then LatPointings++
  ;if(LongPointings MOD 2 eq 1) then LongPointings++
  
  LatStart = -double(LatPointings-1)/2.0*HalfField ;-0.042*deg
  ;LongStart = -double(LongPointings-1)/2.0*HalfField

  cnt=1
  openw, lun,  'print_bulge_survey_grid.dat', /get_lun
  openw, lun2,  'print_bulge_survey_grid.reg', /get_lun
  for o=0L, LongPointings-1 do begin
     for a=0L, LatPointings-1 do begin
        Lat = (LatStart + a*HalfField)/deg
        Lon = (LongMin + o*HalfField)/deg

        EULER, lon, lat, ra, dec, 2
        printf,lun,'00000',cnt,'001',ra,dec,format='(a,i04,a,2f14.6)'
        printf,lun2,'fk5; circle(',ra,',',dec,', ', HalfField/arcmin/60.0,')', format='(a,f14.6,a,f14.6,a,f14.6,a)'
        cnt++
     endfor
  endfor
  close, lun2
  free_lun, lun2
  close, lun
  free_lun, lun
  
end
