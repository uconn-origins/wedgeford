pro b_raytrace,ti,x1,y1,z1,kx,ky,kz,k,totnumph,ctdr,deadphot

common GlobalVars

	;================================================
	;
	;================================================

	rcol=LONG(255*randomu(seed,1))  ; assign a colour to trajectory
	
	smalleststep=1.E-13  ; steps smaller than this are ignored because
		; when kx<0 and cone, we may reintercept current edge,
		; but the calculated distance is slightly too small and 
		; if we insist upon steps>0 then we may forever try to approach
		; this same edge.
		; detect these infinite loops, can't intercept same edge >twice 
		; consequetively.  Use this to set smalleststep, which is otherwise
		; zero, as it should be.		

	large=1.000D250
	tiny=1.000D-250

	itype=0
	numscats=0

 
  ;wset,3
  ;!p.multi=[0,2,2]	

;====================================================
;  itype value corresponds to: 
;====================================================

  itypec=['itype=0','Cyl-1 both -ve',$
    'Cyl-2 +ve and -ve',$
    'Cyl-3 -ve and +ve',$
    'Cyl-4 both +ve',$
    'Cyl-5 ray misses',$
    'Cyl-6 tangent',$
    'Plane-7 only case ',$
    'Cone-8 No intercept',$
    'Cone-9 Tangent',$
    'Cone-10 both -ve',$
    'Cone-11 +ve and -ve',$
    'Cone-12 -ve and +ve',$
    'Cone-13 both +ve',$
    'Cone-14',$
    'Cone-15 ray on surface',$
    'Cone-16']

;====================================================
; 
;====================================================

	vev=tr.arr(*,ti)    ; current cell's vertices elements
	vxv=tr.x(*,ti)   ; current cell's vertices x coords
	vzv=tr.z(*,ti)   ; current cell's vertices z coords

;=============================================================
; establish quantities at outset.
;=============================================================

	steps=0
	bound=1         ; =0 means we're inside grid.
	wabs=0.
	abstau=0.
	;scatoabs=(1.-albedo)/albedo
	tau1=0.d
	tau2=0.d
	rayl=0.d
	currentedge=-100  ; currentedge=-1 at outset, allows us to put source anywhere inside cell,
			  ; let it find it's way initially
	newedge=-100
	nkx=kx
	nky=ky
	nkz=kz
	
	u0=[0.,0.,0.]  ; turn off bulk flows.  Make sure we understand +/- effect on x
	
	method=0
	
	turns=0
	cumf=0
	uvec=0

	inwardflag=0
	dummyflag=0

	cellvisit=0
	lypassedflag=0
	
	flagstuck=0  ; Determine if the photon got stuck.

	prevwabs=1.

  freqx=0.  ; var needs to be defined even if not doing Ly-a...

;=============================================================
; First scattering centre
;=============================================================

	rn=randomu(seed,1)
	tautoscat=-1.*alog(1-rn(0))
	;tautoscat=-1.*alog(1-rn(0))+100000. ; DEBUG, ensures no scats

	tauatscat=tautoscat

;=============================================================
; Propagate ray.
;=============================================================

	;print, 'STARTING PHOTON'
	;print, '...pause...'
	;wait, 3
	pausy=''

	ox1=x1
	oy1=y1
	oz1=z1

	kxprev=kx
	kyprev=ky
	kzprev=kz

	tiprev=ti

;==================================================================================
;==================================================================================
;==================================================================================
;
; 			START PHOTON PROPAGATION LOOP
;
;==================================================================================
;==================================================================================
;==================================================================================


 ; ctdr = 0
  x1sv = fltarr(para.maxscats)
  z1sv = fltarr(para.maxscats)

  while (bound eq 1) do begin
    
      steps=steps+1L
      tau1=tau2
        
      ;=========================================
      ; collect the various opacities for this cell.  
      ; First, dust.
      ;=========================================
             
      scad=cell.sca[ti]*cell.epsi[ti]
      absd=cell.abs[ti]*cell.epsi[ti]
      g=cell.g[ti]
      ;print, scad,absd,g
        
      ;stop
        
      ;========================================
      ;  voigt profile
      ;========================================
          
      sigma=0.0
          
      if(flags.lyman) then begin
          freqdop=cell.dop[ti]       ;                                                        
          hfrac=cell.fh[ti]
          a=0.5*const.freql/freqdop  ; Voigt shape parameter = Lorentz width / Doppler width
            ;displacement from line-centre [in Hz]  
            ; in principle we could add bulk motions here by 
            ;Doppler shifting const.freqcen by 1+(k dot v_bulk)/c    
          freqx=(para.frequency-const.freqcen)/freqdop
          ;phi=voigtprofile(1.,freqx,a)  ; normalized basic voigt profile.  
          phi=voigt(a,freqx)  ; normalized basic voigt profile.  
            ;sigma=hfrac*phi*const.f*const.rootpi*const.esq/(const.me*const.c*freqdop)   ; cm2 per H *nucleus* at freqx. Note hfrac
          sigma=hfrac*phi*const.voigtatomic/freqdop   ; cm2 per H *nucleus* at freqx. Note hfrac
          ;print, 'sigma, phi, freqdop,a :',sigma, phi, freqdop,a
      endif
         
      sca=sigma+scad   ; i.e add dust scattering to lyman alpha (lya will usually dominate)
      ;print, 'sca cm2/H=',sca
      ;=========================================
      ; Calculate line absorption as a pure opacity. PER H NUCLEUS
      ;=========================================
        
      ;B_line_absorption, abslines, freqdop, hfrac
      ;abslines=0.
      abslines=0.
        
      ;=========================================
      ; Calculate opacity due to water PER H NUCLEUS
      ;=========================================
        
      ;abswater=waterfraccell(ti)*watercrossxatwav
      absmol=0.
        
      ;=========================================
      ;  Total abs opacity per H NUCLEUS
      ;=========================================
        
      abs=absd+abslines+absmol
        
      ;=========================================
      ;  
      ;=========================================
      ;print, comx(ti),comz(ti)
      ;print, 'ABS dust, d2g       :',absdust,d2g
      ;print, 'ABS water,waterfrac :',abswater,waterfraccell(ti),h2fraccell(ti)
      ;print, 'SCA Lyman per Hnuc  :',sigma
      ;hisx(steps)=x1
        
      ;  vev=trang(*,ti)   ; current cell's vertices elements
      ;  vxv=trangx(*,ti)   ; current cell's vertices elements
      ;  vzv=trangz(*,ti)   ; current cell's vertices elements
      vev=tr.arr(*,ti)     ; current cell's vertices elements
      vxv=tr.x(*,ti)       ; current cell's vertices elements
      vzv=tr.z(*,ti)       ; current cell's vertices elements
        
      ;=========================================
      ; Check point is inside triangle
      ;=========================================
        	
      ; B_inside_triangle, x1,y1,z1, ox1,oy1,oz1,vxv,vzv
        
      ox1=x1
      oz1=z1
        
      ;----------------------------------
        
      uvec=double([large,large,large]) ; Size three, for times when currentedge=-1
      inttype=double([-1.,-1.,-1.]) ; Size three, for times when currentedge=-1
        
  
      for edg1=0,2 do begin 
    	   ;if(edg1 ne currentedge) then begin  ; intercepts for noncurrent edge.  Saving of 33%
           ;if(edg1 ne currentedge*!pi) then begin  ; intercepts for noncurrent edge.  Saving of 33%
    		  edg2=(edg1+1) mod 3  ; These are actually corners, but that's okay.
    		  ;print, ''
          ; print, 'Current edge :',currentedge
          ;print, 'Current photon position :',x1,z1
          ;print, 'Current photon direction', kx,ky,kz
          ;print, ''
          ; print, 'Cell :',ti,' Edge from vertices :',edg1,edg2
    
    	  x3=double(vxv(edg1))      ; form vector of current cell edge
          x4=double(vxv(edg2))
          z3=double(vzv(edg1))
          z4=double(vzv(edg2))
    
    	  ;==============================================
          ;find intersection of trajectory and each edge
          ;using local.wasp.uwa.edu.au/~pbourke/geometry/lineline2d/
          ; and www.geometrytools.com
          ;==============================================
          b_intersection,x1,z1,kx,ky,kz,x3,z3,x4,z4,ua,large,tiny,itype,dummyflag  ; replace call with code: it is faster. 
          uvec(edg1)=ua
          inttype(edg1)=itype
          ;print, 'itype',itype
          ;print, itypec(itype), itype,ua, ti,k ;,fup
    	;endif 
      endfor
  
      ;==============================================
      ;Find minimum distance
      ;==============================================
  
      ;print, 'DISTANCES :',uvec
  	
      ;print, 'intercepts:',inttype
  
      ;whgtz=where(uvec gt 0. and finite(uvec) eq 1,numwhpos)
      ;whgtz=where(uvec gt smalleststep and finite(uvec) eq 1,numwhpos)
      whgtz=where(uvec gt smalleststep and uvec lt large eq 1,numwhpos)

      if(numwhpos eq 0) then begin
            
         ;	print,'************ Problem, no edge found'
         ; 	print, '    UVEC    :',uvec
         ; 	print, '    INTTYPE :',inttype
         ;       print, '    position:',x1,y1,z1
         ;       print, '    direction:',kx,ky,kz
         ;       print, '    vxv:',vxv
         ;       print, '    vzv:',vzv
         ;       print, '    x3,x4:',x3,x4
         ;       print, '    z3,z4:',z3,z4
 ;         	print, '    Escape, start new photon'
          ; stop

        ; If less than ideal option, try closer cell: 
          whgtz=where(uvec gt 1e-18 and uvec lt large eq 1,numwhpos)
     endif


     if (numwhpos eq 0) then begin
           diag.fup=diag.fup+1L
          bound=0
          goto, escaperoute
      endif

    ; if (z1 lt 0.5) then begin
    ; print, '***** No Problem'
    ;            print, '    UVEC    :',uvec
    ;            print, '    INTTYPE :',inttype
    ;            print, '    position:',x1,y1,z1
    ;            print, '    direction:',kx,ky,kz
    ;            print, '    vxv:',vxv
    ;            print, '    vzv:',vzv
    ;            print, '    x3,x4:',x3,x4
    ;            print, '    z3,z4:',z3,z4
    ;   endif
     ;==============================================
     ;winning edge, cell.  Also find new edge in winning cell
     ;==============================================
      uvecgtz=uvec(whgtz)
      winnergtz=where(uvecgtz eq min(uvecgtz))
     ; print, 'min of uvecgtz',min(uvecgtz),winnergtz,whgtz
      winner=whgtz(winnergtz)  ; winning edge 0,1, or 2
     ; print, 'Of 0, 1 or 2, closest edge is :',winner,' in current cell.'
      winningcell=tr.neigh[winner,ti]
     ; print, 'Current cell :',ti
     ; print, 'Next cell is :',winningcell
      if(n_elements(winningcell) ne 1) then begin
          		print, '***** HORROR!'
          		print, 'FUP distances :',uvec
          		print, 'FUP intercepts:',inttype
          		print, 'FUP directions:',kx,ky,kz
          		print, 'FUP positions:',x1,y1,z1
           print, 'itype',itype
          ;print, ''
          diag.fup=diag.fup+1L
          bound=0
          ;xyouts, x1,z1,'HERE'
          ;if(plotFlags.overplotmain) then oplot, [x1],[z1],psym=1, color=254
          ;stop
          goto, escaperoute
      endif


 
      if(winningcell ne -1) then begin
  		    nedgevec=tr.neigh(*,winningcell)
  		    newedge=where(nedgevec eq ti)
          ;print, 'Current edge in current cell       currentedge :',currentedge, trneigh(*,ti)
          ;print, 'Winning edge in current cell       winner      :',winner
          ;print, 'Moving to new edge in winning cell     newedge :',newedge, trneigh(*,winningcell)
          ;print, ''
          currentedge=newedge
      endif

    ;read, booby

						;================================================
						;  Plot.  different styles
;						;================================================
;
;						if(plotFlags.overplotmain eq 2) then begin
;
;							oplot, [tr.x[0,winningcell],tr.x[1,winningcell]],$
;							[tr.z[0,winningcell],tr.z[1,winningcell]], thick=2
;
;							oplot, [tr.x[1,winningcell],tr.x[2,winningcell]],$
;							[tr.z[1,winningcell],tr.z[2,winningcell]], thick=2
;
;							oplot, [tr.x[2,winningcell],tr.x[0,winningcell]],$
;							[tr.z[2,winningcell],tr.z[0,winningcell]], thick=2
;
;						endif 
;
;						if(plotFlags.overplotmain eq 3) then begin
;
;						;wset, 3
;						;!p.multi=[0,2,1]
;							plot, [tr.x(0,tiprev),tr.x(1,tiprev)],$
;							[tr.z(0,tiprev),tr.z(1,tiprev)], thick=2, $
;							yrange=[min(tr.z(*,tiprev))-0.1*(max(tr.z(*,tiprev))-min(tr.z(*,tiprev))),max(tr.z(*,tiprev))+0.1*(max(tr.z(*,tiprev))-min(tr.z(*,tiprev)))],$
;							xrange=[min(tr.x(*,tiprev))-0.1*(max(tr.x(*,tiprev))-min(tr.x(*,tiprev))),max(tr.x(*,tiprev))+0.1*(max(tr.x(*,tiprev))-min(tr.x(*,tiprev)))], /xs,/ys
;
;							oplot, [tr.x(1,tiprev),tr.x(2,tiprev)],$
;							[tr.z(1,tiprev),tr.z(2,tiprev)], thick=2
;
;							oplot, [tr.x(2,tiprev),tr.x(0,tiprev)],$
;							[tr.z(2,tiprev),tr.z(0,tiprev)], thick=2
;
;							;oplot, [ox1], [oz1], psym=1, color=200
;							oplot, [x1], [z1], psym=2, color=200
;						
;							plot, [tr.x(0,ti),tr.x(1,ti)],$
;							[tr.z(0,ti),tr.z(1,ti)], thick=2, $
;							yrange=[min(tr.z(*,ti))-0.1*(max(tr.z(*,ti))-min(tr.z(*,ti))),max(tr.z(*,ti))+0.1*(max(tr.z(*,ti))-min(tr.z(*,ti)))],$
;							xrange=[min(tr.x(*,ti))-0.1*(max(tr.x(*,ti))-min(tr.x(*,ti))),max(tr.x(*,ti))+0.1*(max(tr.x(*,ti))-min(tr.x(*,ti)))], /xs,/ys
;
;							oplot, [tr.x(1,ti),tr.x(2,ti)],$
;							[tr.z(1,ti),tr.z(2,ti)], thick=2
;
;							oplot, [tr.x(2,ti),tr.x(0,ti)],$
;							[tr.z(2,ti),tr.z(0,ti)], thick=2
;
;							;oplot, [ox1], [oz1], psym=1, color=200
;							oplot, [x1], [z1], psym=2, color=200
;						endif
        
      ;==============================================
      ; Crossing cell, evaluate what this entails -> scatter?
      ;==============================================
        
      xfactor=1.0D  ; These can used to 'put' photon on =0 boundaries etc.
      yfactor=1.0D  ; 
      zfactor=1.0D  ; 
        
      ustep=0.
      ustep=uvec(winner(0))
        
      rayl=rayl+ustep
      dtausca=cell.nh[ti]*sca*grid.xphyscm*ustep    ; eventually use more general expression, rather than density
      dtauabs=cell.nh[ti]*abs*grid.xphyscm*ustep	; d2g already in abs, sca, since sca include ly-alpha
      tau2=tau1+dtausca
   
      ;print, tau1,tau2, dtausca, ':',tauatscat
       
      x1prev=x1   ; capture current before updated.  In the next round, if we error,
      y1prev=y1   ; these variables will capture the previous point to error
      z1prev=z1   
      kxprev=kx
      kyprev=ky
      kzprev=kz    
      tiprev=ti

      if(tau2 gt tauatscat) then begin   ; will we overstep scat point?  These taus are cumulative.
        	;=====================================================
        	;  Scatter.  Test for which of dust/Ly-lpha occurs
        	;=====================================================
                
        	;print, 'scattering'
        	;print, 'Scatter  tau1,tau2,dtau: tauatscat',tau1,tau2,dtau, tauatscat
        	;stop
        	;beep
        
        	;  using tau2 regress x,y to scatter point.  
        
        	ox1=x1   ; use of ox etc perhaps is unnecessary
        	oy1=y1
        	oz1=z1
        
        	; Proper step to scattering centre
        	
        	fracstep=(tauatscat-tau1)/dtausca	; step to scat centre as fraction within cell.
        
        	;   Reorganise this, fewer ops.
        
        	wabs=exp(-1.*abstau)
        
        	ph.j(ti)=ph.j(ti)+wabs*ustep*fracstep      ;*fracstep*(1.-exp(-1.*dtau*fracstep))/(dtau*fracstep) 
        	ph.l(ti)=ph.l(ti)+ustep*fracstep
        
        	ph.fx(ti)=ph.fx(ti)+wabs*ustep*fracstep*kx 
        	ph.fz(ti)=ph.fz(ti)+wabs*ustep*fracstep*kz 
        
        	cellvisit=1.  ; if we scatter this means we must've visited cell.
        			; this var is *only* to limit calls to histogram.  
        			;Nothing to do with ncell.
        	
        	abstau=abstau+dtauabs*fracstep            ; incrememnt by *dust* abs through *total* sca only.
        
        	x1=ox1+kx*ustep*fracstep    ; new position, at scat centre
        	y1=oy1+ky*ustep*fracstep
        	z1=oz1+kz*ustep*fracstep
        
        	currentedge=-100 ; all edges viable, even one we just left.

        	;================================================
        	;  Decide type of scattering.	
        	;================================================
        	
        	r=randomu(seed)
        	dustrat=scad/sca  ; prob of it being dust scattering.
        	dustrat=1.
        	
        	if(r le dustrat) then begin
        	
          	;================================================
          	;  Dust scatter
          	;================================================
        		;print, 'Dust scattering, prob:',dustrat, r
        		
        		b_scatter_dust,g,kx,ky,kz,nkx,nky,nkz,numscats,xrayflag
        		
        	endif else begin
        	
          	;================================================
          	;  Line scatter.  No freq redist.  
          	;  Simple Lya = treat like isotropic dust (g=0.01).
          	;  Or just do phi = acos (2p-1), theta=2pi*p
          	;================================================
             ;b_scatter_dust,0.01,kx,ky,kz,nkx,nky,nkz,numscats
             
             harvest=randomu(seed,2)
 
             kz=(2.0*harvest[0]-1.0)
             
             theta=2*!pi*harvest[1]
             kx=cos(theta)
             ky=Sin(theta)     	
             
        	endelse
        
        	;====================================
        	; correct tau2 and new scat center
        	;====================================
        
        	newti=ti
        	tau2=tau1+dtausca*fracstep
        
        	;print, 'Having stepped, current pos  =  scat point:',tau2,tauatscat
        	
        	;===================================================
        	;  New scattering centre
        	;===================================================
        	
        	rn=randomu(seed,1)
        	tautoscat=-1.*alog(1-rn(0))
        	tauatscat=tau2+tautoscat    ; new scat point
        	tau1=tau2                   ; current position updated to scat centre
        
        	numscats++
        ;print, 'current pos, new scat centre :',tau2,tauatscat
      endif else begin   

        	;===================================================
        	; ELSE (i.e no scatter) move to winning edge ...
        	;===================================================
        	
        	newti=winningcell[0]  ;The anticipated new cell
        
        	nkx=kx  ; default: stay on course
        	nky=ky  ; xceptions are given below
        	nkz=kz
        
        	if(newti eq -1) then begin   
        		
        		; print, 'BOUNDARY DETECTED.  '
        		; 2=top of domain
        		; 3=Outer edge of domain
        		; 4=midplane
        
        		bound=0  ; default.  
        
        		if(tr.edgetype(winner,ti) eq 5) then begin     ; midplane
        			;print,'Reflecting at mid-plane'
        			currentedge=winner
        			nkz=-1.*nkz   ; reflection
        			bound=1       ; still in play
        			newti=ti      ; remain in cell
        			nfreqx=freqx  ; same cell, same doppler.  For systematic flow, need modify by v.dot.k term
        			zfactor=0.    ; Puts us exactly on midplane, over-writes fact that we have reversed nkz and so won't naturally get to midplane.
        		endif
        		
        		if(tr.edgetype(winner,ti) eq 4) then begin     ; r=0.  Meaningful for planar.  Let us reflect.
        			;print,'Reflecting at r=0'
        			currentedge=winner
        			nkx=-1.*nkx   ; reflection
        			bound=1       ; still in play
        			newti=ti      ; remain in cell	
        			nfreqx=freqx  ; same cell, same doppler.  For systematic flow, need modify by v.dot.k term
        			xfactor=0.    ; necessary?
        		endif
        		
        		if(tr.edgetype(winner,ti) eq 2) then begin     ; top  edge only
        			;print, 'exiting at top (or corner)'
        			if(tr.edgetype(winner,ti) eq 2 or tr.edgetype(winner,ti) eq 3) then begin     ; both top and outer edges
        				bound=0
        			endif
        		endif
          
          endif else begin  ; end boundary
          
        		;==============================================================
        		;Put here events contingent upon going to new NON-DOMAIN EDGE cell
        		; i.e that require newti NE -1
        		;==============================================================
          
        		ph.n[newti]++   ; advancement into NEW cell increment.
          
        		;==============================================================
        		;  New freqx due to change in doppler freq of new cell
        		;==============================================================
        		;print, 'old freqx :',freqx
        		;nfreqx=freqx*cell.dop[ti]/cell.dop[newti]  ; +shift due to systematic vel?	
        		;print, 'new freqx :',freqx

          endelse    ; end no-scatter step
          
          ;print, ''
          ;print, 'Crossing cell without scattering: tau1,tau2',tau1,tau2

          wabs=exp(-1.*abstau)             ; this is dust abs weighting 'wabs' as packet enters cell, not yet updated.
	
          ph.j(ti)=ph.j(ti)+wabs*ustep   ;CORRECTION FACTOR -> *(1.-exp(-1.*dtau))/dtau  ?
		                         ;distance weighted intensity contribution exp term is ave. exp^-tau over dtau
          ph.l(ti)=ph.l(ti)+ustep

          ph.fx(ti)=ph.fx(ti)+wabs*ustep*kx 
          ph.fz(ti)=ph.fz(ti)+wabs*ustep*kz 

          cellvisit=0   ; 

          abstau=abstau+dtauabs            ; incrememnt by *dust* abs through *total* sca only.
          ;print, 'No scat. dtauabs, abstau',dtauabs, abstau
          tau1=tau2

          ox1=x1
          oy1=y1
          oz1=z1

          x1=(ox1+kx*ustep)*xfactor   ; the new location
          y1=(oy1+ky*ustep)*yfactor 
          z1=(oz1+kz*ustep)*zfactor 

      endelse

      ;=====================================
      ;  Update quantities, regardless of event.
      ;=====================================

      ;if(sampleel(ti) gt 0 and cellvisit eq 0) then begin

      ;freqxunit=freqx*dopcell(ti)/dopunit

      ;bin=max(where(histbins lt freqxunit, whcnt))
      ;if(whcnt ne 0) then samplehist(bin,sampleel(ti)-1)=samplehist(bin,sampleel(ti)-1)+wabs
      ;print, 'Photon added to cell :',sampleel(ti)-1

      ;cellvisit=1

      ;endif


      ;=====================================
      ;  UPDATE
      ;=====================================

      ti=newti

      kx=nkx
      ky=nky
      kz=nkz

      nkx=kx
      nky=ky
      nkz=kz

      newti=ti
      ;  Transform from 3D to 2D
      ;=====================================

      b_transform, x1,y1,x1trans,y1trans,kx,ky,kxtrans,kytrans

      if (finite(x1trans) eq 0) then begin
 ;         stop
 ;         print,'Transformed to infinity, setting to large.',x1trans,x1
         ctdr = ctdr + 1
         x1trans=large
      endif
  
      x1=x1trans
      y1=y1trans
      kx=kxtrans
      ky=kytrans

      ;print, 'immediately OUT of TRANSFORM kx:',kx

;=====================================
;  Plot trajectory segments
;=====================================

      if(plotFlags.nodevice EQ 0 AND plotFlags.overplotmain EQ 1 AND plotFlags.global EQ 1) then begin
          ;wset,1
          ;print, 'overpotting window'
          oplot, [ox1,x1]*grid.xphys,[oz1,z1]*grid.xphys, color=k*254./totnumph, thick=2.
          ;oplot, [x1]*grid.xphys,[z1]*grid.xphys, psym=8, color=0, symsize=1
          ;if(newti NE -1) then oplot, [cell.comx[newti]]*grid.xphys,[cell.comz[newti]]*grid.xphys,psym=2, color=255
      endif

      escaperoute:
          inwardflag=0
          ; print, ''

     ;===================================
     ;  Is photon packet below attenuation limit for dust absorption?
     ;===================================

;     if(wabs LT para.minwabs) then begin
     if(wabs LE para.minwabs) then begin
 ;         print,'Dim photon (dust+water) prevwabs,wabs,dtauabs :',[prevwabs,wabs, dtauabs]   ;,k,numscats]
          bound=0
          deadphot = deadphot + 1
          ;stop
     endif
     
     if(numscats gt para.maxscats) then begin
          bound=0
          ;print, 'too many scatterings, stopping photon...'	
     endif
     
     if(x1 gt 1d200 or y1 gt 1d200 or z1 gt 1d200) then begin
          bound=0
          deadphot = deadphot + 1
     endif    
     
     if(x1 eq ox1 and y1 eq oy1 and z1 eq oz1) then begin
          if (flagstuck eq 1) then begin
            bound=0
            print, 'Photon stuck twice.'
            deadphot = deadphot + 1
          endif
          flagstuck = 1
     endif     
     ;x1sv[numscat] = x1
     ;z1sv[ctdr] = z1
     ;ctdr = ctdr + 1
     
     prevwabs=wabs
     if steps gt 10000 then begin
      deadphot = deadphot + 1
      bound=0
      print, 'Steps: ',steps,x1,y1,z1
     endif
     
 endwhile

;print, 'End, k, numsctas, wabs, log(wabs) :',[k, numscats, wabs, -alog(wabs)]

end
