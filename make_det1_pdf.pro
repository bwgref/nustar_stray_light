;; Makes DET1 image by projecting PDF of each photon (PDF is selected based on grade).

pro make_det1_pdf,root=root,obsid=obsid,ab=ab,silent=silent, grade_max=grade_max, key=key, livetime=livetime, evtpath=evtpath

  npix=64
  det1w=360
  half=180


  if(n_elements(evtpath) eq 0) then evtpath='event_cl'
  if(n_elements(key) eq 0) then key='key'
  if(n_elements(grade_max) eq 0) then grade_max=4
  if(n_elements(silent) eq 0) then silent=0
  if(n_elements(root) eq 0) then begin
     print,"Give me root directory! root='path to data'" 
     return
  endif
  if(n_elements(obsid) eq 0) then begin
     print,"Give me obsid! e.g. obsid='40010001002'"
     return
  endif

  if(getenv('CALDB') eq '') then begin
     print,'*** ERROR ***'
     print,'Hey, where is NuSTAR $CALDB variable?'
     stop
  endif
  caldb=getenv('CALDB')+'/'
  dirinstrmap=caldb+'data/nustar/fpm/bcf/instrmap/'
  dirpixpos=caldb+'data/nustar/fpm/bcf/pixpos/'
  dirbadpix=caldb+'data/nustar/fpm/bcf/badpix/'
  file=dirpixpos+'nu'+ab+'pixpos20100101v005.fits'

;  pp1=mrdfits(file,1,hh,/silent)
;  pp2=mrdfits(file,2,hh,/silent)
;  pp3=mrdfits(file,3,hh,/silent)
;  pp4=mrdfits(file,4,hh,/silent)

  pp=[[mrdfits(file,1,hh,/silent)], $
      [mrdfits(file,2,hh,/silent)], $
      [mrdfits(file,3,hh,/silent)], $
      [mrdfits(file,4,hh,/silent)]]

  ;print,n_elements(pp[*,0])
  ;help, pp[*,0],/str
  ;stop

  fpd = lonarr(npix/2.,npix/2.,4) 
  mask_fpd = lonarr(npix/2.,npix/2.,4)+1 
  mask_fpd[0,*,*]=0
  mask_fpd[31,*,*]=0
  mask_fpd[*,0,*]=0
  mask_fpd[*,31,*]=0

  det1=DBLARR(det1w,det1w)
  evt = root+'/'+obsid+'/'+evtpath+'/nu'+obsid+ab+'01_cl.evt'
  print,'Read '+evt

  if(~file_test(evt)) then begin
     print, 'Not found:', evt
     stop
  endif
  print,'Read ',evt

  copy = root+'/'+obsid+'/'+evtpath+'/nu'+obsid+ab+'01_cl_pdf.evt'
  file_copy, evt, copy, /overwrite
  a = mrdfits(evt, 1, h)
  naxis2 = SXPAR(h,'NAXIS2',COUNT=cnt)
  livetime = SXPAR(h,'LIVETIME',COUNT=cnt)

  ; modify header to add new column
  var=0.0
  fxbaddcol, index, h, var, 'weight', 'pixel PDF weight'

  ; modify data to add new column
  wfield = {weight: 0.0}
  addfield = replicate(wfield,naxis2)
  a = struct_addtags(a,addfield)
  
  ; read badpix info 
  for id=0,3 do begin
     table=readfits(evt,header,EXTEN_NO=(3+id),/SILENT)
     badrawx=TBGET(header, table, 'RAWX', /NOSCALE)
     badrawy=TBGET(header, table, 'RAWY', /NOSCALE)
     for r=0,n_elements(badrawx)-1 do mask_fpd[badrawx[r],badrawy[r],id]=0
  endfor

  if(1) then begin
     addrow=0L
     for i=0L,50 do begin
     ;for i=0L,N_ELEMENTS(a)-1 do begin
     if (silent ne 1) then $
        print,String(13b),(i*1.0/(N_ELEMENTS(a)-1))*100.0,format='(A,f12.4,"%",$)'
     ix=a[i].rawx
     iy=a[i].rawy
     grade=a[i].grade
     if not (mask_fpd[ix,iy,a[i].det_id] eq 1) then continue

     count=0
     det_id=a[i].det_id

     ;; find pixpos configuration for this event (RAWX/Y and grade)
     ii=where(pp[*,det_id].rawx eq ix and pp[*,det_id].rawy eq iy and $
              pp[*,det_id].grade eq grade and pp[*,det_id].ref_det1x ne -1,count)

     if(count eq 1) then begin
        if (finite(total(pp[ii,det_id].pdf))) then begin
           det1[pp[ii,det_id].ref_det1x:pp[ii,det_id].ref_det1x+6,$
                pp[ii,det_id].ref_det1y:pp[ii,det_id].ref_det1y+6] += pp[ii,det_id].pdf
           for py=0L,6 do for px=0L, 6 do begin
              index=px+py*7
              if (pp[ii,det_id].pdf[index] gt 0.0) then begin
                 b=a[i]
                 b.det1x=pp[ii,det_id].ref_det1x+px
                 b.det1y=pp[ii,det_id].ref_det1y+py
                 b.weight=pp[ii,det_id].pdf[index]
                 a=struct_append(a,b)
                 addrow++
              endif
           endfor
        endif
     endif

  ;;    case a[i].det_id of
  ;;       0: begin
  ;;          ii=where(pp1.rawx eq ix and pp1.rawy eq iy and pp1.grade eq grade and pp1.ref_det1x ne -1,count)  
  ;;          if(count eq 1) then begin
  ;;             if (finite(total(pp1[ii].pdf))) then begin
  ;;                det1[pp1[ii].ref_det1x:pp1[ii].ref_det1x+6,$
  ;;                     pp1[ii].ref_det1y:pp1[ii].ref_det1y+6] += pp1[ii].pdf
  ;;                for py=0L,6 do for px=0L, 6 do begin
  ;;                   index=px+py*7
  ;;                   if (pp1[ii].pdf[index] gt 0.0) then begin
  ;;                      b=a[i]
  ;;                      b.det1x=pp1[ii].ref_det1x+px
  ;;                      b.det1y=pp1[ii].ref_det1y+py
  ;;                      b.weight=pp1[ii].pdf[index]
  ;;                      a=struct_append(a,b)
  ;;                      addrow++
  ;;                   endif
  ;;                endfor
  ;;             endif
  ;;          endif
  ;;       end
  ;;       1: begin
  ;;          ii=where(pp2.rawx eq ix and pp2.rawy eq iy and pp2.grade eq grade and pp2.ref_det1x ne -1,count)      
  ;;          if(count eq 1) then begin
  ;;             if (finite(total(pp2[ii].pdf))) then begin
  ;;                det1[pp2[ii].ref_det1x:pp2[ii].ref_det1x+6,$
  ;;                     pp2[ii].ref_det1y:pp2[ii].ref_det1y+6] += pp2[ii].pdf
  ;;                for py=0L,6 do for px=0L, 6 do begin
  ;;                   index=px+py*7
  ;;                   if (pp2[ii].pdf[index] gt 0.0) then begin
  ;;                      b=a[i]
  ;;                      b.det1x=pp2[ii].ref_det1x+px
  ;;                      b.det1y=pp2[ii].ref_det1y+py
  ;;                      b.weight=pp2[ii].pdf[index]
  ;;                      a=struct_append(a,b)
  ;;                      addrow++
  ;;                   endif
  ;;                endfor
  ;;             endif
  ;;          endif
  ;;       end
  ;;       2: begin
  ;;          ii=where(pp3.rawx eq ix and pp3.rawy eq iy and pp3.grade eq grade and pp3.ref_det1x ne -1,count)      
  ;;          if(count eq 1) then begin
  ;;             if (finite(total(pp3[ii].pdf))) then begin
  ;;                det1[pp3[ii].ref_det1x:pp3[ii].ref_det1x+6,$
  ;;                     pp3[ii].ref_det1y:pp3[ii].ref_det1y+6] += pp3[ii].pdf
  ;;                for py=0L,6 do for px=0L, 6 do begin
  ;;                   index=px+py*7
  ;;                   if (pp3[ii].pdf[index] gt 0.0) then begin
  ;;                      b=a[i]
  ;;                      b.det1x=pp3[ii].ref_det1x+px
  ;;                      b.det1y=pp3[ii].ref_det1y+py
  ;;                      b.weight=pp3[ii].pdf[index]
  ;;                      a=struct_append(a,b)
  ;;                      addrow++
  ;;                   endif
  ;;                endfor
  ;;             endif
  ;;          endif
  ;;       end
  ;;       3: begin
  ;;          ii=where(pp4.rawx eq ix and pp4.rawy eq iy and pp4.grade eq grade and pp4.ref_det1x ne -1,count)   
  ;;          if(count eq 1) then begin
  ;;             if (finite(total(pp4[ii].pdf))) then begin
  ;;                det1[pp4[ii].ref_det1x:pp4[ii].ref_det1x+6,$
  ;;                     pp4[ii].ref_det1y:pp4[ii].ref_det1y+6] += pp4[ii].pdf
  ;;                for py=0L,6 do for px=0L, 6 do begin
  ;;                   index=px+py*7
  ;;                   if (pp4[ii].pdf[index] gt 0.0) then begin
  ;;                      b=a[i]
  ;;                      b.det1x=pp4[ii].ref_det1x+px
  ;;                      b.det1y=pp4[ii].ref_det1y+py
  ;;                      b.weight=pp4[ii].pdf[index]
  ;;                      a=struct_append(a,b)
  ;;                      addrow++
  ;;                   endif
  ;;                endfor
  ;;             endif
  ;;          endif
  ;;       end
  ;;       else: PRINT, 'Not one through four'
  ;;    endcase 
  endfor
  endif

  print,naxis2,addrow
  sxaddpar,h,'NAXIS2',addrow
  ii=where(a.weight gt 0.0,cnt)
  if (cnt eq 0) then stop
  modfits,copy,a[ii],h,exten_no=1

  mkhdr, hdr_det1_long, 3, [det1w,det1w]
  SXADDPAR, hdr_det1_long, 'LIVETIME', livetime, ' Exposure, s.'
  SXADDPAR, hdr_det1_long, 'OBSID', obsid, ''
  SXADDPAR, hdr_det1_long, 'MODULE', ab, ' Module'
        
  writefits,obsid+ab+'.'+key+'.det1.fits',det1[*,*],hdr_det1_long
     
end
