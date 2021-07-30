
pro b_photon_initial_stellar,kx,ky,kz,k,x1,y1,z1,totnumph,ti

common GlobalVars

;======================================================================
;
;  k is photon counter.
;
;======================================================================

;=======================================================================
;  given chosen position, find nearest cell com and use that instead.
;=======================================================================

  chpos=[0.005,0.005]  ; roughly Sun radius.  IN AU !!! not code units
  
  diff=(cell.comx*grid.xphys-chpos(0))^2.+(cell.comz*grid.xphys-chpos(1))^2.
  whdiff=where(diff eq min(diff))
  fixti=whdiff(0)

  ti=fixti  ; initial cell

  y1=0.0D
  x1=0.01*cell.comx[fixti]   ; place us in centre of cell to start.
  z1=0.01*cell.comz[fixti] 

;========================
;  Angles
;========================

if(flags.randomizeangle eq 1) then begin

	theta=acos(randomu(seeda1,/double)*(cos(const.degtorad*para.theta2)-cos(const.degtorad*para.theta1))+cos(const.degtorad*para.theta1))  

endif else begin

	theta=acos((1.*(k+1.)/totnumph)*(cos(const.degtorad*para.theta2)-cos(const.degtorad*para.theta1))+cos(const.degtorad*para.theta1))  

endelse

;	angle.phi=0.
;
;	kx=double(sin(angle.theta)*cos(angle.phi))  ; may want to redefine theta, phi?
;	ky=double(sin(angle.theta)*sin(angle.phi))
;	kz=double(cos(angle.theta))

  kx=double(sin(theta))
  ky=0.
  kz=double(cos(theta))

;print, 'Stellar kx,ky,kz,k,x1,y1,z1,totnumph,ti:',kx,ky,kz,k,x1,y1,z1,totnumph,ti


end
