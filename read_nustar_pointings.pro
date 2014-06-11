; read NuSTAR pointings in format: "RA Dec"
pro read_nustar_pointings,file, pnt_ra, pnt_dec,pnt_lon=pnt_lon,pnt_lat=pnt_lat,obsid=obsid

READCOL,file,F='A,F,F',obsid,pnt_ra,pnt_dec,delimiter=" ",count=count,nlines=file_nlines
print,'read ',count,' lines out of ',file_nlines,' from ',file
if(count eq 0) then begin 
   print,'No lines read from ',file
   stop
endif
EULER, pnt_ra, pnt_dec, pnt_lon, pnt_lat, 1
end
