pro read_source_catalog,src_file,src_name,src_ra,src_dec,src_flag,index=index, nsrc=nsrc

READCOL,src_file,F='A,F,F,B',src_name0,src_ra0,src_dec0,src_flag0, delimiter=",",count=src_count,nlines=src_file_nlines
print,'read ',src_count,' lines out of ',src_file_nlines,' from ',src_file
if(src_count eq 0) then begin 
   print,'No lines read from ',src_file
   return
endif
index=where(src_flag0 eq 1,selected)
nsrc=selected

case size(src_name, /n_dimen) of
  0: src_name = src_name0
  1: src_name = [src_name, src_name0]
endcase

case size(src_ra, /n_dimen) of
  0: src_ra = src_ra0
  1: src_ra = [src_ra, src_ra0]
endcase

case size(src_dec, /n_dimen) of
  0: src_dec = src_dec0
  1: src_dec = [src_dec, src_dec0]
endcase

case size(src_flag, /n_dimen) of
  0: src_flag = src_flag0
  1: src_flag = [src_flag, src_flag0]
endcase


if(selected eq 0) then begin 
   print,'No source selected '
   return
endif else begin
   print,'Selected ',selected
endelse
end
