;;
;; xx varies in nthe range 0-63
;;
function det2mm, xx

  hgap = 0.15                   ; detector half gap [mm]  

  return, (xx+1-32.5)*0.6048 + hgap*sign(xx+1-32.5)

end
