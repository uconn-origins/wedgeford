pro main2
 
common GlobalVars
 
 pi=!dpi


	;plots, cell.comx(fixti)*grid.xphys,cell.comz(fixti)*grid.xphys, psym=2, color=255, thick=1, symsize=1  ; plot star

;==========================================
;  Loop over photons starts here
;==========================================

  ctdr = 0L
  deadphot = 0L
	freqx=0.
	width=100.

	print, 'Starting rad trans'
	fup=0L
  imgcnt=0L

	compfrac=0.1
	compn=1.*LONG(compfrac*para.totnumph)

for k=0L,LONG(para.totnumph-1) do begin

	if(k/compn eq LONG(k/compn)) then begin
	
		print, 'Photon :',k+1,' of ',para.totnumph
		; print, 'Infphot: ',ctdr
		; print, 'Deadphot: ',deadphot
		
	endif

	if(para.isrfg eq 0.0) then begin
	;=====================================================================
	;  Stellar field: set photon start direction 
	;  Stellar field: set photon start location *within* start cell..
	;=====================================================================
	
	     b_photon_initial_stellar,kx,ky,kz,k,x1,y1,z1, para.totnumph,ti ; potentially random angle
         
	endif else begin
  
	;=====================================================================
	;  Interstellar field: set photon start direction 
	;  Interstellar field: set photon start location *within* start cell..
	;=====================================================================

	     b_photon_initial_ISRF,kx,ky,kz,k,x1,y1,z1, para.totnumph,ti ; potentially random angle
  
	endelse

	;==============================================
	;  Propagate photon/ray.
	;==============================================

	b_raytrace,ti,x1,y1,z1,kx,ky,kz,k,para.totnumph,ctdr,deadphot

	;==============================================
	;  screen dump?
	;==============================================

;  img=tvrd(true=1)
;  acceptStr=''
  
 ; filly='img3_'+STRTRIM(STRING(imgcnt, FORMAT='(I04)'),2)
 ; read, acceptStr, PROMPT='Get image?'
 ; if(acceptStr='y') then begin
  ;  imgcnt++
    ;write_png, 'Output/'+filly+'.png', img
 ;   write_jpeg, 'Output/'+filly+'.jpg', img, true=1
 ;
  ;endif


;skiprt:

endfor

print, 'RT done.',diag.fup,' messed up photons. Out of', para.totnumph

end
