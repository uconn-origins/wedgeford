pro B_transform, x1,y1,x1rot,y1rot,kx,ky,kxrot,kyrot

;=====================================================
;  Given a line segment of photon trajectory between cell faces,
; we transform kx,ky,kz into new local coord system
;=====================================================

	;print, 'ORIG: x,y:',x1,y1
	;print, 'ORIG: kx,ky,kx^2.+ky^2. :',kx,ky, kx^2.+ky^2.
       
	magr=sqrt(x1*x1+y1*y1)

	if(magr GT 0.) then begin
		stheta=double(y1/magr)
		ctheta=double(x1/magr)
	endif else begin       ; i.e if we're at r=0 then do pi rotation.
		print, '********* TRANSFORM: in r=0 statement'
		ctheta=-1.0D
		stheta=0.0D
	endelse

	;print, 'Theta, sin, cos :',acos(ctheta),stheta,ctheta

	;x1rot=ctheta*x1+stheta*y1  ;  Actually, xrot is just magr, and yrot=0.
	;y1rot=ctheta*y1-stheta*x1

	x1rot=magr
	y1rot=0.0D

	kxrot=double(ctheta*kx+stheta*ky)
	kyrot=double(ctheta*ky-stheta*kx)

	;print, 'ROTATED: xrot,yrot:',x1rot,y1rot
	;print, 'ROTATED: kxrot,kyrot, sq:',kxrot,kyrot, kxrot^2.+kyrot^2.
	;print, ''

end
