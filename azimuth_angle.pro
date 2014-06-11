; Calculate Azimuth anglep from North to East                                                                                                    

FUNCTION AZIMUTH_ANGLE, ra, dec, ra0, dec0

  vec_src=radec2xyz([ra, dec]) ; Vector to a bright source nearby                                                                                

  vec_pnt=radec2xyz([ra0, dec0]) ; Pointing vector                                                                                               

  vec_pole = [0, 0, 1] ; North pole                                                                                                              

  vec_pole_proj = vec_pole-(total(vec_pole*vec_pnt))*vec_pnt ; North pole projected onto a plane normal to the pointing vector                    

  vec_src_proj = vec_src-(total(vec_src*vec_pnt))*vec_pnt ; Bright source vector projected onto a plane normal to the pointing vector             

  acos_value = acos(total(vec_pole_proj*vec_src_proj)/(sqrt(total(vec_pole_proj*vec_pole_proj))*sqrt(total(vec_src_proj*vec_src_proj))))

  if ra gt ra0 then begin 

     return, acos_value 

  endif else begin 
     
     return, 2.*!PI-acos_value 
     
  endelse

end
