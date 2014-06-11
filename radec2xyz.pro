; Convert RA, DEC to XYZ                                                                                                                         

function radec2xyz, radec

  ra = double(radec[0])
  dec = double(radec[1])
  xyz = [cos(dec)*cos(ra), cos(dec)*sin(ra), sin(dec)]
  return, xyz

end

