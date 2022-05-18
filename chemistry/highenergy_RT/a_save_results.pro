pro a_save_results

common GlobalVars

      date=systime(0)
      dow=strmid(date,0,3)+'_'
      mon=strmid(date,4,3)+'_'
      day=strtrim(strmid(date,8,2)+'_',1)
      time=strmid(date,11,5)+'_'
      year=strmid(date,20,4)
      
      dateStr=mon+day+year
      print, dateStr

;====================================
;  Save result recast on original disk coords
;====================================
 
     fillysave0='./Output/'+strtrim(file.save,2)+'_'+strtrim(file.disk,2)+'_'+dateStr+'.txt'
     print, 'Saving recast results in :', fillysave0
 
openw, 4, fillysave0
  
  printf, 4, 'Results recast onto the original points (before any grid was added)'
  printf, 4, ''
  printf, 4, 'column 1 : x [AU]'
  printf, 4, 'column 2 : z [AU]'
  printf, 4, 'column 3 : n_ph [cm^-3]'
  printf, 4, ''
  for i=0L,LONG(n_elements(orig.xp)-1) do begin
    if(orig.xp[i] LE grid.xupper AND orig.zp[i] LE grid.zupper) then begin
      printf,4,orig.xp[i]*grid.xphys,orig.zp[i]*grid.xphys,orig.nph[i]
    endif else begin
    printf,4,orig.xp[i]*grid.xphys,orig.zp[i]*grid.xphys,0
    endelse
  endfor

close, 4

;====================================
;  Brute force save of all IDL variables.
;====================================
      
      fillysave1='./Output/'+strtrim(file.save,2)+'_'+strtrim(file.disk,2)+'_'+dateStr+'_all.sav'
      ;fillysave2='Output/'+strtrim(file.save,2)+'_'+strtrim(file.disk,2)+'_'+datestr+'_vars.sav'
      
      print, 'Saving entire IDL memory :',fillysave1;,fillysave2
    
      save, filename=fillysave1, /all
      ;save, filename=fillysave2, /variables


end