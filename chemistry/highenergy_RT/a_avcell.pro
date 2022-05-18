pro A_avcell, scadustv, absdustv

common GlobalVars
;common trangblock2
;common domainblock2
;common histblock2
;common waterblock2
;common colblock2

	avcell=0.*dcell
	h2colcell=0.*dcell
	h2colinitcell=0.*dcell
	hcolcell=0.*dcell
	hcolinitcell=0.*dcell

	large=1.E35
	tiny=1.E-15
	fup=0L
	steps=0L
	overplotmain=0

;=====================================
; Direction.  Can be anything, 
; i.e vertical, horizontal, towards star.
;=====================================

	kx=0.
	ky=0.
	kz=1.

	kmag=sqrt(kx*kx+ky*ky+kz*kz)

	kx=kx/kmag
	ky=ky/kmag
	kz=kz/kmag

;=====================================
; Loop over cells, but only those for
; whom we sample radiation field.
;=====================================

print, 'A_avcell : Av along dirn: ',[kx,ky,kz]

print, 'DO WE NEED TO INVOLVE dust:gas ratio for this Av? For consistency, yes, surely.'

for i=0L, LONG(ntr-1) do begin

if(sampleel(i) gt 0) then begin

	;print, 'accepting point :',i

	;=====================================
	;  Initial
	;=====================================
	
	x1=comx(i)
	z1=comz(i)
	y1=0.	

	vev=trang(*,i)    ; current cell's vertices elements
	vxv=trangx(*,i)   ; current cell's vertices x coords
	vzv=trangz(*,i)   ; current cell's vertices z coords

	currentedge=-100  ; not one edge
	newedge=-100
	steps=0L
	
	tautot=0.
	h2col=0.
	hcol=0.
	h2colinit=0.
	watercol=0.
	watercolinit=0.
	ohcol=0.
	ohcolinit=0.

	bound=1	
	ti=i

	racol=255.*randomu(seed)
		
	;=====================================
	;  Propagate photon from starting point
	;=====================================

	while (bound eq 1) do begin

		;print,i,bound,steps,tautot

		uvec=double([large,large,large]) 
		inttype=double([-1.,-1.,-1.]) 

		vev=trang(*,ti)   ; current cell's vertices elements
		vxv=trangx(*,ti)   ; current cell's vertices elements
		vzv=trangz(*,ti)   ; current cell's vertices elements

		;=====================================
		;  Loop over edges
		;=====================================

		for edg1=0,2 do begin

			if(edg1 ne currentedge) then begin  ; intercepts for noncurrent edge

				edg2=(edg1+1) mod 3  ; These are actually corners, but that's okay.

				x3=double(vxv(edg1))      ; form vector of current cell edge
				x4=double(vxv(edg2))
				z3=double(vzv(edg1))
				z4=double(vzv(edg2))

				B_intersection,x1,z1,kx,ky,kz,x3,z3,x4,z4,ua,large,tiny,itype,dummyflag  ; replace call with code: it is faster. 

				uvec(edg1)=ua
	
			endif 

		endfor

	;=====================================
	;  find closest edge and associated cell.
	;=====================================

	whgtz=where(uvec gt tiny and finite(uvec) eq 1,numwhpos)

	if(whgtz eq [-1]) then begin
		print,'************ Problem, no edge found'
		bound=0
		goto, escaperoute
	endif

	uvecgtz=uvec(whgtz)              
	winnergtz=where(uvecgtz eq min(uvecgtz))
	winner=whgtz(winnergtz)  ; winning edge 0,1, or 2

	winningcell=trneigh(winner,ti)

	if(n_elements(winningcell) ne 1) then begin	
		print, '************ NO next cell found. At',i,'of',ntr-1,nsample
		fup=fup+1L
		bound=0
		xyouts, x1,z1,'HERE'
		goto, escaperoute
	endif

	;==============================================
	; Crossing cell, increment tau_tot
	;==============================================

	xfactor=1.  ; this is purely for r=0 situation.
	yfactor=1.  ; this is purely for r=0 situation.
	zfactor=1.  ; this is purely for r=0 situation.

	ustep=0.
	ustep=uvec(winner(0))
 
	;==============================================
        ; Increment Av columns etc.  NOTE: Av includes dust to gas ratio
        ; Note that tau is scat+abs extinction.
        ;==============================================
       
	dtau=dcell(ti)*d2gcell(ti)*(scadustv+absdustv)*xphyscm*ustep    ; eventually use more general expression, rather than density
	tautot=tautot+dtau
 
	dh2col=dcell(ti)*h2fraccell(ti)*xphyscm*ustep   
        h2col=h2col+dh2col

	dhcol=dcell(ti)*hfraccell(ti)*xphyscm*ustep   
        hcol=hcol+dhcol

	dwatercol=dcell(ti)*waterfraccell(ti)*xphyscm*ustep 
	watercol=watercol+dwatercol

	dwatercolinit=dcell(ti)*waterfraccellinit(ti)*xphyscm*ustep 
	watercolinit=watercolinit+dwatercolinit
	
	dohcol=dcell(ti)*ohfraccell(ti)*xphyscm*ustep 
	ohcol=ohcol+dohcol

	dohcolinit=dcell(ti)*ohfraccellinit(ti)*xphyscm*ustep 
	ohcolinit=ohcolinit+dohcolinit
  
	;==============================================
        ; Boundary event?
        ;==============================================

	newti=winningcell(0)  ;The anticipated new cell

	nkx=kx  
	nky=ky  
	nkz=kz

	;==============================================
	; Boundary event?
	;==============================================

	if(newti eq -1) then begin
		
		;=======================================
		; print, 'BOUNDARY DETECTED.  '
		; 2=top of domain
		; 3=Outer edge of domain
		; 4=midplane
		;=======================================

		bound=0  ; default.  

		if(edgetype(winner,ti) eq 5) then begin     ; midplane
			print,'Reflecting at mid-plane'
			currentedge=winner
			nkz=-1.*nkz   ; reflection
			bound=1       ; still in play
			newti=ti      ; remain in cell
			zfactor=0.    ; Puts us exactly on midplane, over-writes fact that we have reversed nkz and so won't naturally get to midplane.
		endif
		
		if(edgetype(winner,ti) eq 4) then begin     ; r=0.  Meaningful for planar.  Let us reflect.
			print,'Reflecting at x=0'
			currentedge=winner
			nkx=-1.*nkx   ; reflection
			bound=1       ; still in play
			newti=ti      ; remain in cell	
			xfactor=0.    ; necessary?
		endif
		
		if(edgetype(winner,ti) eq 2) then begin     ; top  edge only
	
			avcell(i)=tautot

			watercolcell(i)=watercol
			ohcolcell(i)=ohcol
			h2colcell(i)=h2col
			hcolcell(i)=hcol
			watercolinitcell(i)=watercolinit
			ohcolinitcell(i)=ohcolinit
			h2colinitcell(i)=h2colinit

			
			;print,'comx,comy',comx(i),comz(i)
			;print,'Tau_ext :',tautot

		endif

	endif 

;=====================================
;  UPDATE
;=====================================

	ox1=x1
	oy1=y1
	oz1=z1

	x1=(ox1+kx*ustep)*xfactor   ; the new location
	y1=(oy1+ky*ustep)*yfactor 
	z1=(oz1+kz*ustep)*zfactor 

	ti=newti

	kx=nkx
	ky=nky
	kz=nkz

	if(overplotmain EQ 1) then begin

		;oplot, [ox1,x1],[oz1,z1], color=racol, thick=1.
		
	endif

	steps++	

	endwhile

		;print, 'Photon exited (bound=0)',i

	escaperoute:

endif
endfor

end
