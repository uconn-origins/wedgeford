
pro b_photon_initial_ISRF,kx,ky,kz,k,x1,y1,z1,totnumph,ti

common GlobalVars

;======================================================================
;
;  Note: k is photon counter, not mod of dirn vector.
;
;======================================================================

;x=randomu(seeda1,/double)*grid.xupper
;y=0.0
;z=grid.zupper*0.9

  chpos=[grid.xupper*sqrt(randomu(seeda1,/double)),grid.zupper*0.9]  ; along upper boundary, in code units.
  ;chpos=[0.000,0.005]  ; roughly Sun radius.  IN AU !!! not code units
  
  diff=(cell.comx-chpos[0])*(cell.comx-chpos[0])+(cell.comz-chpos[1])*(cell.comz-chpos[1])
  whdiff=where(diff eq min(diff))
  fixti=whdiff(0)

  ti=fixti  ; initial cell

  y1=0.0D
  x1=cell.comx[fixti]   ; place us in centre of cell to start.
  z1=cell.comz[fixti] 

;============================
; Angles
;============================

if(flags.randomizeangle eq 1) then begin

  theta=acos(-1.0*randomu(seeda1,/double))

endif else begin

  theta=acos(-1.0*k/totnumph)

endelse

;theta=!pi*0.75

	phi=2.*!pi*randomu(seeda1,/double)

	kx=double(sin(theta)*cos(phi))  ; may want to redefine theta, phi?
	ky=double(sin(theta)*sin(phi))
	kz=double(cos(theta))


;print, 'ISRF: kx,ky,kz,k,x1,y1,z1,totnumph,ti:',kx,ky,kz,k,x1,y1,z1,totnumph,ti
end
