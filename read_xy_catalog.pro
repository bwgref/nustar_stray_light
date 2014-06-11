pro read_xy_catalog,src_file,src_ra,src_dec

READCOL,src_file,F='F,F',src_ra,src_dec,delimiter=" ",count=src_count,nlines=src_file_nlines
print,'read ',src_count,' lines out of ',src_file_nlines,' from ',src_file
if(src_count eq 0) then begin 
   print,'No lines read from ',src_file
   return
endif
end
