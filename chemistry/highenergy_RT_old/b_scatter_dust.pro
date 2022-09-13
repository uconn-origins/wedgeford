pro b_scatter_dust,g,u,v,w,nu,nv,nw,numscats,xrayrun
    
    
    ;============================================
    ; u,v,w are kx,ky,kz dirn vectors
    ;============================================
      
    harvest=randomu(seed,2)
        
    if xrayrun eq 0 then begin
      cosalpha=(0.5/g)*((1.+g*g)-((1.-g*g)/(1.-g+2.*g*harvest(0)))*((1.-g*g)/(1.-g+2.*g*harvest(0))))
    endif else begin 
      thetaang = weighted_choice_b(harvest(0))
      cosalpha=cos(thetaang)
;      print, 'in xrayrun',cosalpha,thetaang,harvest(0)
    
    endelse
    
    sinalpha=sqrt(1-cosalpha*cosalpha)
    phi=!pi*(2.*harvest(1)-1.)
    
    ;=======================================================
    ;  If numscats eq 0 force phi to send photon downwards, towards disk.
    ;=======================================================
    
    ;if(numscats eq 0 and forced eq 1) then 
    
    ;=======================================================
    ;
    ;=======================================================
    
    sinphi=sin(phi)
    cosphi=cos(phi)
    
    wterm=1/(sqrt(1-w*w))    ; The w term that appears below, saves comp time
    
    nu=  u*cosalpha  -  (v*cosphi  +  w*u*sinphi)  *  sinalpha  *  wterm 
    
    nv=  v*cosalpha  + (u*cosphi  -  v*w*sinphi)  *  sinalpha  *  wterm
    
    nw=  w*cosalpha +   sinalpha  *  sinphi  /  wterm
    
    ;print, 'New dirn mod:',sqrt(nu*nu+nv*nv+nw*nw)
    
    ;nu=abs(nu)  ; forces scatter in +ve x dirn, for debugging
end
