
pro make_badpix,dmask_fp1,dmask_fp2

  badpix_path=getenv('NUPLAN_AUXIL')+'/BadPix/20100101v002'

  time = 78913712.

  det2raw,dmask_fp1,dmask_fp1d
  det2raw,dmask_fp2,dmask_fp2d
  
for t=0,3 do begin
   file='data_userbadpix_A'+String(t,Format='(i1)')+'.dat'
   openw, lun,  file, /get_lun
   for i=0,31 do begin
      for j=0,31 do begin
         
         if(dmask_fp1d[i,j,t] gt 0) then begin
            printf, lun, time, $
                    i, j, $
                    format = '(d10.0, 2i6, "   4")'
         endif
      endfor
   endfor
   close, lun
   free_lun, lun
endfor

for t=0,3 do begin
   file='data_userbadpix_B'+String(t,Format='(i1)')+'.dat'
   openw, lun,  file, /get_lun
   for i=0,31 do begin
      for j=0,31 do begin
         
         if(dmask_fp2d[i,j,t] gt 0) then begin
            printf, lun, time, $
                    i, j, $
                    format = '(d10.0, 2i6, "   4")'
         endif
      endfor
   endfor
   close, lun
   free_lun, lun
endfor

; Brian's code - add edges
files_to_adjust = file_search('data_userbadpix_??.dat')

nfiles = n_elements(files_to_adjust)

FOR i = 0, nfiles - 1 DO begin
   openw, lun,  files_to_adjust[i], /get_lun, /append
   FOR j = 0, 31 DO BEGIN
      badx = 31
      bady = 31
      printf, lun, time, $
                badx, j, $
                format = '(d10.0, 2i6, "   4")'
      printf, lun, time, $
                j, bady, $
                format = '(d10.0, 2i6, "   4")'
   ENDFOR
   close, lun
   free_lun, lun
ENDFOR
; end of Brian's code

; copy all needed files to the current directory, and run shell script
aux=['header_userbadpix*.txt','rm_history.txt',$
'header_all.txt', 'primary_fpma.txt',$
'primary_fpmb.txt', 'columns_userbadpix.txt',$
'create_badpix_fpma.txt', 'create_badpix_fpmb.txt']

for i=0,n_elements(aux)-1 do spawn,'cp '+badpix_path+'/'+aux(i)+' .'


ar=['a','b']
for i=0,1 do begin
   runme='create_badpix_fpm'+ar(i)+'.txt'
   openr, lun, runme, /get_lun
   rows = File_Lines(runme)
   cmds=STRARR(rows)
   ReadF, lun, cmds
   close, lun
   free_lun, lun
   spawn,strjoin(cmds, ';'), log
   openw, lun,  'badpix_fp'+ar(i)+'.log', /get_lun
   printf,lun,log
   close, lun
   free_lun, lun
endfor

for i=0,n_elements(aux)-1 do spawn,'rm '+aux(i)
spawn, 'rm data_userbadpix_??.dat'

end
