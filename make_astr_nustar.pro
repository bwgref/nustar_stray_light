; makes astr structure for NuSTAR snapshot
pro make_astr_nustar, center_ra, center_dec, PA, astr=astr, oa=oa

  if(n_elements(oa) eq 0) then oa=[38., 38.] ; "default" position of OA (+3mm,+3mm)
  dr = !PI/180.
  ;print,'make_astr_nustar: Generate rotation matrix for PA=',PA

  rot_mat = [ [ cos(PA*dr), -sin(PA*dr)], $ ;Rotation matrix
              [ sin(PA*dr),  cos(PA*dr)] ] 

  crpix=oa
  cdelt=[12.3/3600.0, 12.3/3600]
  make_astr,astr, CD=rot_mat, DELTA = cdelt, CRPIX = crpix, CRVAL = [center_ra,center_dec]

  ; it forces FK5, but sometimes we need galactic coords
  ;RADECSYS = 'FK5', EQUINOX = 2000.0
end

