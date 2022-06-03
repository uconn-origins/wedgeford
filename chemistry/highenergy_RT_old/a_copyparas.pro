pro a_copyparas,paraFile

common GlobalVars

  date=systime(0)
  ;dow=strmid(date,0,3)+'_'
  mon=strmid(date,4,3)+'-'
  day=strtrim(string(LONG(strmid(date,8,2))),1)+'_'
  timehr=strmid(date,11,2)+'h'
  timemn=strmid(date,14,2)
  year=strmid(date,20,4)+'-'
  
  dateStr=year+mon+day+timehr+timemn
  print, '* Making backup of parameter file.  The date :',dateStr

  parafilly='./History/paras_'+strtrim(file.save,2)+'_'+strtrim(file.disk)+'_'+dateStr+'.txt'

  file_copy, './'+paraFile ,parafilly, /verbose, /overwrite
  
  openw, 2, parafilly, /APPEND
  printf,2,' '
  printf,2,'===================== SUMMARY =========================='
  printf,2,' '
  
  if(flags.useabs EQ 0) then printf,2, 'You suppressed dust absorption entirely!' 
  if(flags.usesca EQ 0) then printf,2, 'You suppressed dust scattering entirely!' 
  if(para.constabs NE 0) then printf,2, 'Enforced const/homog dust abs cm^2/H :',para.constabs
  if(para.constsca NE 0) then printf,2, 'Enforced const/homog dust sca cm^2/H :',para.constsca
  if(para.consteps NE 0) then printf,2, 'Enforced const/homog dust epsilon :',para.constsca
  if(para.constg NE 0) then printf,2, 'Enforced const/homog dust g=<cos> :',para.constsca
  
  
  if(flags.lyman NE 0) then printf,2, 'This was a Lyman alpha run'
  if(flags.castgrid NE 0) then begin
    printf,2, 'You used a cartesian background grid, base Nr x Nz:',grid.nxg, grid.nzg  
    printf,2, '... with ',flags.castnested,' levels of nesting.'
  endif
  if(flags.castradial NE 0) then begin
    printf,2, 'You used a radial power-law r^-p z^-q background grid, base Nr x Nz:',grid.nxg, grid.nzg 
    printf,2, '... with p,q=',grid.px,grid.pz,' power-laws.'
  endif
  
;  if(flags.lyman NE 0) then printf,2, ''
;  if(flags.lyman NE 0) then printf,2, ''
;  if(flags.lyman NE 0) then printf,2, ''
;  if(flags.lyman NE 0) then printf,2, ''
  
  printf,2,' '
  printf,2,'===================== SAMPLE OF INPUT DISK FILE =========================='
  printf,2,' '
  
  openr, 3, 'MyDisk/'+strtrim(file.disk,2)+'.txt'
  lineStr=''
  
  for i=0,20 do begin
  
    readf,3, lineStr
    printf,2,lineStr
  
  endfor
  
  close, 3
  close, 2



end

