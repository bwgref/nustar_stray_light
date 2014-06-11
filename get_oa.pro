DEPRECATED 

pro get_oa,oa,pos=pos
  common nuplan, nu, status, sources, target

   if(nu.oa[0] gt max(nu.xpos_array) or $
      nu.oa[0] lt min(nu.xpos_array) or $
      nu.oa[1] gt max(nu.ypos_array) or $
      nu.oa[1] lt min(nu.ypos_array)) then begin
      print,'Error: OA is out of FOV'
      widget_control, status.mainid, /destroy
      stop
   endif

  xpos_near = Min(Abs(nu.xpos_array - oa[0]), xpos_index)
  ypos_near = Min(Abs(nu.ypos_array - oa[1]), ypos_index)

  dx=1
  dy=1

  if(nu.xpos_array[xpos_index] gt oa[0]) then dx=-1
  if(nu.ypos_array[ypos_index] gt oa[1]) then dy=-1

  pos=fltarr(2)
  pos[0]=(oa[0]-nu.xpos_array[xpos_index+dx])/(nu.xpos_array[xpos_index]-nu.xpos_array[xpos_index+dx])*(-dx)+(xpos_index+dx)
  pos[1]=(oa[1]-nu.ypos_array[ypos_index+dy])/(nu.ypos_array[ypos_index]-nu.ypos_array[ypos_index+dy])*(-dy)+(ypos_index+dy)

end
