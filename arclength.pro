; in radian                                                                                                                                      

FUNCTION arclength, ra1, dec1, ra2, dec2

return, acos(cos(!PI/2.-dec1)*cos(!PI/2.-dec2)+sin(!PI/2.-dec1)*sin(!PI/2.-dec2)*cos(ra1-ra2))

END

