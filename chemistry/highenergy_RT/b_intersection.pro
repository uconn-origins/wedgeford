pro b_intersection,x1,z1,kx,ky,kz,x3,z3,x4,z4,ua,large,tiny,itype,inwardflag

	; x1=current pos, x2=x1+direction vector, such that |r2-r1|=1 
	; x3=edge point one, x4=edge point two, thus making edge vector
	; similarly for y.

	;x3=5.
	;y3=3.
	;x4=6.0
	;y4=8.

	;x1=-9.
	;y1=0.

	;kx=0.8  ;left/right
	;ky=0.2  ; up/down
	;kz=0.0  ; into/out of page

	inwardflag=0

	tsol=large
	itype=-1
	t1=large
	t2=large

	quad=0

	;=========================================
	;  setup for cone, in 3D
	;=========================================

	p=[x1,0D,z1]   ; photon point
	dtemp=[kx,ky,kz]   ;[-1,-0.424,1]                 ; photon direction
	d=dtemp/sqrt(dtemp(0)*dtemp(0)+dtemp(1)*dtemp(1)+dtemp(2)*dtemp(2)) ; normalised, should already be
						; normal when in main program

	;=========================================
	;setup for cylinder, is done in 2D (x,y) 
	;=========================================

	p2d=[x1,0.0D]   ; in 2D, for cylinder routine, z coord not important, and y=0 due to rotation
	dtemp2d=[kx,ky]   ;[-1,-0.424,1]                 ; photon direction
	d2d=dtemp2d/sqrt(dtemp(0)*dtemp(0)+dtemp(1)*dtemp(1)+dtemp(2)*dtemp(2)) ; 2D vector normalised to full 3D version
	; this ensures that solution for cylinder is compatible with 3D dirn vector.


	a=[0.0D,0.0D,1.0D]   ; cone axis, never changes
	r4=[x4,0.0D,z4]
	r3=[x3,0.0D,z3]    

	gradient=(z4-z3)/(x4-x3)
         
	intercept=z4-gradient*x4

	v=[0.0D,0.0D,intercept]   ; vertex of cone

	dif=r4-r3
	gamma=(z4-z3)/sqrt(total(dif*dif))  ;internal angle of cone.

	; x3=x4 test for cylinder is simpler, and lets us put the above in he
	; its respective routine.  Savings.  

;================================================================
; Check to see if it's a cylinder, run cylinder routine
;================================================================

;if (gamma eq 1.) then begin
if (x3 eq x4) then begin

	ddel=d2d(0)*p2d(0)+d2d(1)*p2d(1)
	dsq=d2d(0)*d2d(0)+d2d(1)*d2d(1)
	delsq=p2d(0)*p2d(0)+p2d(1)*p2d(1)

	cdelta=ddel*ddel-dsq*(delsq-x3*x3)

	t1=0.
	t2=0.

	; print, 'CYLINDER CASE. radius',x3
	;=====================================
	;  plot some vectors
	;=====================================

	;loadct, 39
	;plot, [p(0),p(0)+2.*d(0)],[p(2),p(2)+2.*d(2)],$
	; xrange=[-20,20],/xs, yrange=[-20,20],/ys, /iso, thick=5
	;oplot, [x3,x4],[-1000.,1000.], color=200
	;oplot, -1.*[x3,x4],[-1000.,1000.], color=200
	;xyouts, p(0),p(2),'P',alignment=0.5, charsize=2.

	;=====================================

	if(cdelta gt 0.) then begin

		;print,'		Cyl. ray hits cylinder at two points'

		;itype=0	
		t1=(-1.0D*ddel-sqrt(cdelta))/dsq
		t2=(-1.0D*ddel+sqrt(cdelta))/dsq

		quad=1

		if(t1 LE 0.0D and t2 LE 0.0D) then begin

			itype=1
			tsol=large
			;print, '	  Cyl. entirely behind photon', t1,t2
			goto, out
		endif

		if(t1*t2 lt 0.0D) then begin

			;print, '    Cyl. One distance +ve , one -ve'

			if(t1 gt t2 ) then begin
				itype=2
				tsol=t1
				goto, out
			endif else begin 
				itype=3
				tsol=t2
				goto, out
			endelse

		endif else begin  ; end -/+, begin +/+
		
			;print, '	  Cyl. Both +ve'
			itype=4
			tsol=min([t1,t2])
			goto, out
		endelse

	endif else begin

	if(cdelta lt 0.0D) then begin

		itype=5
		tsol=large
		;print, '  	Cyl. ray misses cylinder'
		goto, out

	endif

	if(cdelta eq 0.0D) then begin
		itype=6
		tsol=-1.0D*ddel/dsq
		;print, '	  Cyl. ray tangent to cylinder'
		goto, out
	endif

	endelse

	print, ' ****** Cyl. failed to characterise. cdelta :',cdelta
	;stop
	   tsol=large
	
	goto, out

endif  ;================= end cylinder, otherwise do plane

;=================================================
; check for Plane case.
;=================================================

if(z3 eq z4) then begin    ; plane

	; print, 'PLANE CASE. horizontal at z=',z3

	itype=7
	tsol=(z3-z1)/kz
	;print, '    Plane distance, virgin:',tsol

	; if(tsol le 0.0D or kz eq 0.0D) then tsol=large
	if(kz eq 0.0D) then tsol=large

	; print, '    Plane distance, tested :',tsol
	goto, out

endif

;=====================================
;  plot some vectors for cone case
;=====================================

	;loadct, 39
	;plot, [p(0),p(0)+2.*d(0)],[p(2),p(2)+2.*d(2)],$
	;	 xrange=[-20,20],/xs, yrange=[-20,20],/ys, /iso, thick=5
	;oplot, [v(0),v(0)+100.*sqrt(1.-gamma*gamma)],[v(2),v(2)+100.*gamma], color=200
	;oplot, [v(0),v(0)-100.*sqrt(1.-gamma*gamma)],[v(2),v(2)+100.*gamma], color=200
	;oplot, [v(0),v(0)+100.*sqrt(1.-gamma*gamma)],[v(2),v(2)-100.*gamma], color=200
	;oplot, [v(0),v(0)-100.*sqrt(1.-gamma*gamma)],[v(2),v(2)-100.*gamma], color=200
	;xyouts, v(0),v(2),'V',alignment=0.5, charsize=2.
	;xyouts, p(0),p(2),'P',alignment=0.5, charsize=2.

	;=====================================
	;
	;=====================================

  m=a#a-gamma*gamma*[[1.0D,0.0D,0.0D],[0.0D,1.0D,0.0D],[0.0D,0.0D,1.0D]]

  delt= p-v
; ;print, 'Cone. delt=p-v :',delt
; ;print, 'cone. matrix m:',m
; ;print, 'cone. photon dirn d:',d
;
  c0a= delt # m # delt  ;del m del  ; matrix multiply
  c1a=    d # m # delt  ;d m del
  c2a=    d # m # d     ;d m d

  ;print, 'Cone. Coefficients c2,c1,c0',c2a(0),c1a(0),c0a(0)
  ;print, 'Cone. x3,x4 coords',x3,x4
  ;print, 'Cone. z3,z4 coords',z3,z4

; c0=1.5401139e-14 ;c0a(0)
; c1=1.8701113e-07 ; c1a(0)
; c2=0.66463459 ; :c2a(0)

  c0= c0a(0)
  c1= c1a(0)
  c2= c2a(0)
	;=====================================
	; the following cases are ordered by probability, so if sequence
	; which be exited as soon as possible, on average.
	;=====================================
	
        t1=large
	t2=large

	term=c1*c1-c0*c2

;	print, 'Term :',term

   if(abs(c2) gt tiny) then begin   ; Implies Quadratic, c2 NE 0, but either neg or pos.

	if(term lt 0.0D) then begin
		;print, 'Cone. NO INTERCEPT'
		itype=8
		tsol=large
		goto, out
	endif

	if(term eq 0.0D) then begin
	;	print, 'Cone. TANGENT'
		itype=9
		tsol=-1.*c1/c2
		nsol=1
		goto, out
	endif

	;print, 'Quadratic'
         quad=2
	t1=(-1.0D*c1+sqrt(term))/c2
	t2=(-1.0D*c1-sqrt(term))/c2  ; t1<=t2  ??
 ;       t1=(-1.0D*c1+2.0*sqrt(term))/(2.0*c2)
 ;       t2=(-1.0D*c1-2.0*sqrt(term))/(2.0*c2)  ; LIC trying this based on eq in Appendix B
	
	if(t1 LE 0.0D and t2 LE 0.0D) then begin ; start -/-
		;print, 'Cone. No intercept.'
		itype=10
		tsol=large
		goto, out
	endif

          ;print,'t1: ',t1,' t2: ',t2
          ;print,'C1: ',c1,' C2: ',c2,' C0: ',c0
          ;print,'kz: ',kz
          ;print,'whatever:',(z1+kz*t1-intercept)*gradient
	if(t1*t2 lt 0.0D) then begin 

	;	print, 'Cone. One distance +ve , one -ve'

		if(t1 gt t2 ) then begin
			itype=11
			tsol=t1
			if((z1+kz*tsol-intercept)*gradient le 0.0D) then tsol=large
			goto, out
		endif else begin 
			itype=12
			tsol=t2
			if((z1+kz*tsol-intercept)*gradient le 0.0D) then tsol=large
			goto, out
		endelse
	endif else begin  ; end -/+, begin +/+
	;	print, 'Cone. Both +ve'
		;if(total(finite([t1,t2])) lt 2.0D) then begin
		 ; print, t1,t2
		;endif
		itype=13
		tsol=min([t1,t2])
               ; if t1 eq 0.0 then tsol = t2
                ;if t2 eq 0.0 then tsol = t1
		if((z1+kz*tsol-intercept)*gradient le 0.0D) then tsol=max([t1,t2])
		goto, out
	endelse

	;print, 'Cone.  given that there should be TWO solns, we failed to find them'
	stop
	
endif 

if(abs(c2) lt tiny) then begin   ; Linear eqn, on solution.

	if(abs(c1) gt tiny and abs(c0) gt tiny) then begin	
		itype=14
		tsol=-0.5D*(c0/c1)
		goto, out
	endif

	if(abs(c1) lt tiny and abs(c0) lt tiny) then begin
		;print, 'Ray lies on cone, and therefore we have no intercept.'
		itype=15
		tsol=large
		goto, out
	endif

	if(abs(c0) lt tiny and abs(c1) lt tiny) then begin
		itype=16
		tsol=large
		goto, out
	endif
 
	;print, 'Cone.  Given that there should be ONE soln, we failed to find it.'
	stop
	
endif
	
;==========================
;
;==========================

out:

	;==========================
	;
	;==========================
	; do we need to distinguish intercepts between
	; original cone and image cone?  In general no, because another cone
	;  will win before the image cone does.
	;  But, perhaps an issue does arise.  From cone vector we can 
	; identify image/original since Z coord changes side of vertex 
	; when we cross to image cone


	;==========================
	;  plots
	;==========================

	;ps1=p+d*t1(0)
	;ps2=p+d*t2(0)
	;ps3=p+d*tsol;*[tsol,tsol,tsol]
	;print, 'distances :',t1,t2
	;print, 'intersection at ;',ps1
	;print, 'intersection at ;',ps2
	;plots,ps1[0],ps1[2], psym=2
	;plots,ps2[0],ps2[2], psym=2
	;plots,ps3[0],ps3[2], psym=2, symsize=3

	;
	; 
	;

	

;	if(inwardflag eq 1) then begin
;
;		if(quad gt 0) then begin
;
;			tsol=max([t1,t2]) ; assumes the supposed zero is always a +ve roudnign error
;			if(t1*t2 LT 0.0D AND max([t1,t2]) LT tiny) then tsol=large
;			if(t1 lt 0.0D and t2 lt 0.0D) then tsol=large
;	
;		endif

		;if(quad eq 0) then begin

		;endif
		;print,'* Inward intercept: t1, t2, tsol, quad :', t1,t2, tsol,quad

;	endif

ua=tsol
inwardflag=0
 ;prin0t, 'ua=',ua,inwardflag
    ;    print, '-----'

end
