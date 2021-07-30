pro a_compute_packet_luminosity

common GlobalVars

;==========================================
;  L_star=1.  Evrything in photon number density
;  Note upper hemisphere of star =0.5Lstar
;  solid angle 2pi.  
;==========================================

	solidfrac=0.5*(cos(const.degtorad*para.theta1)-cos(const.degtorad*para.theta2))  ; only one hemisphere is done.

	ph.packetlumin=solidfrac/para.totnumph    ; packet luminosity as *fraction* of stellar luminosity as given in parameter file
	
	print, '* Frac solid angle coverage :',solidfrac
	print, '* Photon packet luminosity  (Lstar=1):',ph.packetlumin

end
