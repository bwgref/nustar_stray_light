
; e.g.
;root='/Data/NuSTAR/data/mini-survey'
;obsid_array=['40010001002', '40010002001', '40010003001', '40010004001', '40010005001', '40010006001']

pro make_badpix_region,root=root,obsid=obsid,silent=silent, ds9reg=ds9reg, exclude=exclude, firstrun=firstrun, regA=ds9regA, regB=ds9regB, evtsdir=evtsdir
  det1w=360
  hotpix_max=1000
  e1=2.
  e2=10.
  cha1=(e1-1.6)/0.04
  cha2=(e2-1.6)/0.04
  mod_name=['A','B']
  npix=64
  hotpix_raw=lonarr(npix,npix,2)

  print,'make_badpix_region'

  if(n_elements(firstrun) eq 0) then firstrun=0 else firstrun=1
  if(n_elements(exclude) eq 0) then property=1 else property=0
  if(n_elements(ds9reg) eq 0) then begin
     print,"Give me ds9 polygon region file in degrees! ds9reg='path to polygon.reg'" 
     return
  endif
  
  if(n_elements(silent) eq 0) then silent=0
  if(n_elements(border) eq 0) then border=10
  if(n_elements(root) eq 0) then begin
     print,"Give me root directory! root='path to data'" 
     return
  endif
  if(n_elements(obsid) eq 0) then begin
     print,"Give me obsid! e.g. obsid='40010001002'"
     return
  endif


  det1=LONARR(det1w,det1w,2)
  files_evt = root+'/'+obsid+['/'+evtsdir+'/nu'+obsid+'A01_cl.evt','/'+evtsdir+'/nu'+obsid+'B01_cl.evt']
  for t=0,1 do begin
     print,'Read ',files_evt[t]
     table=readfits(files_evt[t],header,EXTEN_NO=1,/SILENT)
     grade=TBGET( header, table, 'GRADE', /NOSCALE)
     chan=TBGET( header, table, 'PI', /NOSCALE)
     det_id=TBGET( header, table, 'DET_ID', /NOSCALE)
     det1x=TBGET( header, table, 'DET1X', /NOSCALE)-1
     det1y=TBGET( header, table, 'DET1Y', /NOSCALE)-1
     xsky=TBGET( header, table, 'X', /NOSCALE)
     ysky=TBGET( header, table, 'Y', /NOSCALE)
     rawx=TBGET( header, table, 'RAWX', /NOSCALE)
     rawy=TBGET( header, table, 'RAWY', /NOSCALE)
     nx = SXPAR(Header,'TLMAX13',COUNT=A)
     ny = SXPAR(Header,'TLMAX14',COUNT=A)
     crval1 = SXPAR(Header,'TCRVL13',COUNT=A)
     crval2 = SXPAR(Header,'TCRVL14',COUNT=A)
     crpix1 = SXPAR(Header,'TCRPX13',COUNT=A)
     crpix2 = SXPAR(Header,'TCRPX14',COUNT=A)
     cdelt1 = SXPAR(Header,'TCDLT13',COUNT=A)
     cdelt2 = SXPAR(Header,'TCDLT14',COUNT=A)
     ctype1 = SXPAR(Header,'TCTYP13',COUNT=A)
     ctype2 = SXPAR(Header,'TCTYP14',COUNT=A)

     make_astr,astr, DELTA = [cdelt1, cdelt2], $
               CRPIX = [crpix1,crpix2], CRVAL = [crval1, crval2], CTYPE=[ctype1,ctype2]
               ;RADECSYS = 'FK5', EQUINOX = 2000.0
     mkhdr, hdr, 3, [nx,ny]
     putast, hdr, astr, CD_TYPE=2

     filereg=ds9reg
     if(t eq 0) then begin
        if not (n_elements(ds9regA) eq 0) then filereg=ds9regA        
     end else begin
        if not (n_elements(ds9regB) eq 0) then filereg=ds9regB        
     endelse
     print,filereg
     imsize=[ny,ny]
     ima=lonarr(nx,ny)
     openr, ds9, filereg, /get_lun
     string=''
     WHILE NOT EOF(ds9) DO BEGIN ; main loop
        readf, ds9, string
        IF strmid(string, 0, 8) NE 'polygon(' THEN goto, skip   
; extract coordinate components
        coors = strmid(string,8,strlen(string)-9)
        coors = float(strsplit(coors,',',/extract))
        npoint = n_elements(coors)/2
        ind = indgen(npoint)
        px = coors[ind*2]
        py = coors[ind*2+1]
        merge_vector, px, px[0] ; to close the loop
        merge_vector, py, py[0]
skip:
     ENDWHILE                   ; end of main loop
     free_lun, ds9
     
     ad2xy,px,py,astr,pdx,pdy
     
     fpd = lonarr(npix/2.,npix/2.,4) 
     total=0L
     for i=0L,N_ELEMENTS(xsky)-1 do begin
        xx=xsky(i)-1
        yy=ysky(i)-1
        
        if(InsidePolygon (xx, yy, pdx, pdy) eq property) then begin
           fpd(rawx(i),rawy(i),det_id(i))++
           total++
        endif else begin
           if(xx ge 0 and yy ge 0 and xx lt nx and yy lt ny) then begin
              ima(xx,yy)++
           endif
        endelse
        
     endfor

     if(silent ne 1) then begin 
        print,'inside polygon: fp'+mod_name[t],total
        writefits,'masked_SKY_reg.fp'+mod_name[t]+'.img',ima,hdr
     endif

  raw2det,fpd,fp
  hotpix_raw[*,*,t]=fp
  writefits,'masked_RAW_reg.fp'+mod_name[t]+'.img',fp
  

  endfor

   ; check is there badpix-region matrix
   filesav='region_badpix.sav'
   if(file_test(filesav) and firstrun eq 0) then begin
      print,'make_bappix_region: RESTORE ',filesav
      hotpix_raw_save=hotpix_raw
      RESTORE, filesav
      hotpix_raw[*,*,0]+=hotpix_raw_save[*,*,0]
      hotpix_raw[*,*,1]+=hotpix_raw_save[*,*,1]
   endif else begin
      print,'make_bappix_region: SAVE ',filesav
   endelse
   save, filename=filesav, hotpix_raw, root, obsid
end
