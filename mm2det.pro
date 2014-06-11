function mm2det, mm

  hgap = 0.15                   ; detector half gap [mm]  

  return, (mm-hgap*sign(mm))/0.6048 -1 +32.5

end
