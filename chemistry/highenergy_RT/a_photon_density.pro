pro a_photon_density

common GlobalVars

;==========================================================
;   jcellnorm.
; From jcell is distance weighted Weight (wabs), code units. 
; lcell is integrated pathlength of all ncell visiting packets
; jcellnorm is total (normalised) luminosity entering cell (volume) per second.
; fxcellnorm is x flux as fraction of jcellnorm
;==========================================================

	ph.jnorm=0.*cell.nh
	ph.nphnorm=0.*cell.nh
	ph.fxnorm=0.*cell.nh
	ph.fznorm=0.*cell.nh

	wlgtz=where(ph.n gt 0.,numwh)
	if(numwh gt 0) then begin

		;ph.jnorm(wlgtz)=para.starlumin*ph.j(wlgtz)		
;	  ph.jnorm[wlgtz]=para.starlumin*(ph.j[wlgtz]/ph.l[wlgtz])  
    ph.jnorm[wlgtz]=ph.j[wlgtz];/ph.l[wlgtz]  
		ph.fxnorm[wlgtz]=ph.fx[wlgtz]/ph.j[wlgtz]                 ; unitless. kx weighted average.  Other factors cancel.
		ph.fznorm[wlgtz]=ph.fz[wlgtz]/ph.j[wlgtz]

	endif

;stop

	wnan=finite(ph.jnorm)
	whnan=where(wnan eq 0.,numnotfin)

	if(numnotfin gt 0) then ph.jnorm(whnan)=0.

;=================================================
;  Photon density in cgs
;=================================================

	print, '+ Max jcellnorm :',max(ph.jnorm)
	
	;ph.nph=(para.starlumin/const.c)*ph.jnorm/(grid.xphyscm*grid.xphyscm*tr.vol)
  ;ph.nph=(para.starlumin*ph.j/const.c)/(tr.vol*grid.xphyscm^2.)
 ;ph.nph=(para.starlumin*ph.j/ph.l/const.c)/(tr.vol*grid.xphyscm^2.)
	ph.nph=(para.starlumin*ph.packetlumin/const.c)*ph.jnorm/(tr.vol*grid.xphyscm^2.)
		
;		print, ph.nph
;		print, ph.jnorm
;		print, ph.j
;    print, tr.vol
;    print, tr.vol*(grid.xphyscm*grid.xphyscm)/para.starlumin
    
;=================================================
; Photon flux ph cm-2 s-1
;=================================================

	ph.fx=const.c*ph.nph*ph.fxnorm
	ph.fz=const.c*ph.nph*ph.fznorm

; check this line later LIC (should it be ph.fx and ph.fz?)
; changed second ph.fx to ph.fz

  ph.a=(sqrt(ph.fx^2.+ph.fz^2.))/(ph.nph*const.c)

;=================================================
; 
;=================================================
	
	;print, lylum/(4.*!pi*2.25E26),' Lyman ph cm^-2 s^-1  ... at 1AU'
	;print, lylum/(4.*c*!pi*2.25E26),' Lyman ph cm^-3  ... at 1AU'
	
	;endif
;stop

end
