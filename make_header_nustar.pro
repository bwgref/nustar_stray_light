; makes header for NuSTAR snapshot
; bitpix: 1 - Byte, 2 - 16 bit integer, 4 - float, 3 - Long
pro make_header_nustar, center_ra, center_dec, PA, astr=astr, hdr=hdr, bitpix=bitpix, oa=oa
  make_astr_nustar, center_ra, center_dec, PA, astr=astr, oa=oa
  mkhdr, hdr, bitpix, [64,64]  
  putast, hdr, astr, CD_TYPE=2
end
