pro a_casting_main, xp, zp, nh, epsi, sca, abss, g, fh, tg, fwater, foh

common GlobalVars

;common trangblock2
;common flagblock2

;====================================================
;  Initial casting of user-defined points,
;  then unstructured PP points, and finally a 
;  variety of grids
;  
;  the values returned are the values on the nodes.
;====================================================

if(strtrim(file.disk,2) EQ 'none') then begin

  print, 'Creating generator points ...'
  
  ;a_invent_object, xpf, zpf, nhf
  
  ;plot, xpf,zpf, psym=3

endif else begin

	print, 'Getting generator points ...'

  x=dindgen(1)
  z=dindgen(1)
  
  dotrianginterp=1  ; background grids interpolated using triangulate, instead of nearest original data point
  
;======================================================
;  Disk file
;======================================================

  print, 'Opening :','./MyDisk/'+strtrim(file.disk,2)+'.txt'
  ;openf, 2, './MyDisk/'+strtrim(file.disk,2) 
  filly='MyDisk/'+strtrim(file.disk,2)+'.txt'


  readcol,filly,index,xpf,zpf,nhf,epsif,scaf,absf,gf,fhf,tgf,fwaterf,fohf, /SILENT  

   ; print,index[0],xpf[0],zpf[0],nhf[0],epsif[0],scaf[0],absf[0],gf[0],fhf[0],tgf[0],fwaterf[0],fohf[0]
  
  print, 'Number of points in ', strtrim(file.disk,2),' =',max(index)+1

endelse


;=== recast into code units. 
  
  xpfc=xpf/grid.xphys
  zpfc=zpf/grid.xphys
  
  
  orig={xp:xpfc,zp:zpfc,$  ; in code units
        nh:nhf,$
        nph:nhf*0.,$
        fx:nhf*0.,$
        fy:nhf*0.,$
        fz:nhf*0.}   ; need to record these original points for recasting onto new grid at the very end.

;==== only keep those that fall in our domain
  
  whgood=where(xpfc LE grid.xupper AND zpfc LE grid.zupper,numgood)
  
  print, 'Of these ', numgood, ' fall in specified domain.'
  
  xpl=xpfc[whgood]
  zpl=zpfc[whgood]
  
  nhl=nhf[whgood]
  epsil=epsif[whgood]
  scal=scaf[whgood]
  absl=absf[whgood]
  gl=gf[whgood]
  fhl=fhf[whgood]
  tgl=tgf[whgood]
  fwaterl=fwaterf[whgood]
  fohl=fohf[whgood]
 
 ;============= note inner gap is in code units.
 
  grid.xgap=min(xpl)
 
 print, 'Implied gap at [AU]:',grid.xgap*grid.xphys
 
  
 
 
 
 
;=====================  excise hgigh altitude points

;ang=zpl/xpl
;;
;whgood=where(ang LE 0.6,numgood)
;
;  xpl=xpl[whgood]
;  ypl=ypl[whgood]
;  zpl=zpl[whgood]
;  nhl=nhl[whgood]
;
  nh=nhl
  xp=xpl
  zp=zpl

  epsi=epsil
  sca=scal
  abss=absl
  g=gl
  fh=fhl
  tg=tgl
  fwater=fwaterl
  foh=fohl

;=====================  PP

;if(grid.npoisson gt 0) then begin

;	print, 'Casting Poisson PP points.  Loop for more points.' 

;endif

;===========================================================
;  Cast a fixed grid, so cells have maximum size, helps plotting. 
;===========================================================

	if(flags.castgrid eq 1) then begin

  print, 'Casting (possibly nested) grid ...'
  print, 'flags: castgrid, num levels', flags.castgrid, flags.castnested
    
		if(flags.xf EQ 0 and flags.zf EQ 0) then begin
		  print, 'Any nested grids aligned to low-left corner'
		  a_castgrid_nested, xg, zg
		  ;stop
    endif else begin
      print, 'Any nested grids are centered to specified focus point:', flags.xf, flags.zf
      a_castgrid_nested_focus, xg, zg
      ;stop
 		endelse
 		
 	  methodlabel = 'linear'
 	  maxinvdist = 6
		;=======================================
		;assign densities etc. to these points
    ;=======================================
     
    print, 'Assigning densities to grid points...'
  
    nhg=0.*xg
 
      if(dotrianginterp EQ 1) then begin
      
       triangulate, xp, zp, tr
       
       zervec_xp = findgen(100)/101.0
       zervec_zp = fltarr(100)
       zervec_nh = dblarr(100)

       for jz = 0,99 do begin
         junk = min(sqrt((xp-zervec_xp(jz))^2+(zp-zervec_zp(jz))^2),minind)
         zervec_nh(jz) = nh(minind)
         ;print,xp(minind)*grid.xphys,zp(minind)*grid.xphys,nh(minind)
       endfor
 
       zervec_xp_bk = fltarr(100)+0.95
       zervec_zp_bk = findgen(100)/121.0
       zervec_nh_bk = dblarr(100)
       
       for jz = 0,99 do begin
         junk = min(sqrt((xp-zervec_xp_bk(jz))^2+(zp-zervec_zp_bk(jz))^2),minind)
         zervec_nh_bk(jz) = nh(minind)
         ;print,xp(minind)*grid.xphys,zp(minind)*grid.xphys,nh(minind)
       endfor 

       zervec_xp_bk = fltarr(100)+1.05
       ; print,transpose([zervec_xp_bk,zervec_zp_bk,zervec_nh_bk])
       zervec_zp = fltarr(100)-0.05
       nhgII = [nh,zervec_nh ,zervec_nh_bk]
       xpII = [xp,zervec_xp,zervec_xp_bk]
       zpII = [zp,zervec_zp,zervec_zp_bk]
       
       ;print,[transpose(zervec_xp),transpose(zervec_zp),transpose(zervec_nh)]
       
       triangulate, xpII, zpII, trII
       
       nhg=GRIDDATA(xpII, zpII, nhgII, xout=xg, yout=zg, method=methodlabel, triangles=trII)
       ;hg=GRIDDATA(xp, zp, nh, xout=xg, yout=zg, method=methodlabel, triangles=tr)
       epsig=GRIDDATA(xp, zp, epsi, xout=xg, yout=zg, method=methodlabel, triangles=tr)
       scag=GRIDDATA(xp, zp, sca, xout=xg, yout=zg, method=methodlabel, triangles=tr)
       absg=GRIDDATA(xp, zp, abss, xout=xg, yout=zg, method=methodlabel, triangles=tr)
       gg=GRIDDATA(xp, zp, g, xout=xg, yout=zg, method=methodlabel, triangles=tr)
       fhg=GRIDDATA(xp, zp, fh, xout=xg, yout=zg, method=methodlabel, triangles=tr)
       tgg=GRIDDATA(xp, zp, tg, xout=xg, yout=zg, method=methodlabel, triangles=tr)
       fwaterg=GRIDDATA(xp, zp, fwater, xout=xg, yout=zg, method=methodlabel, triangles=tr)
       fohg=GRIDDATA(xp, zp, foh, xout=xg, yout=zg, method=methodlabel, triangles=tr)
       
       
       ; Check if density values are negative or finite.
       badvals_nh = where(finite(nhg) eq 0.0)
       negvals_nh = where(nhg lt 0.0)
       
       ;print,badvals_nh
       ;print,negvals_nh 
       ;if badvals_nh ne -1 then begin
       ;  print,'w'
       ;endif
       
       methodzerzg = 'nearestneighbor'
       if negvals_nh(0) ne -1 then begin
         nhg(negvals_nh)=GRIDDATA(xpII, zpII, nhgII, xout=xg(negvals_nh), yout=zg(negvals_nh), method=methodzerzg, triangles=trII)
         epsig(negvals_nh)=GRIDDATA(xp, zp, epsi, xout=xg(negvals_nh), yout=zg(negvals_nh), method=methodzerzg, triangles=tr)
         scag(negvals_nh)=GRIDDATA(xp, zp, sca, xout=xg(negvals_nh), yout=zg(negvals_nh), method=methodzerzg, triangles=tr)
         absg(negvals_nh)=GRIDDATA(xp, zp, abss, xout=xg(negvals_nh), yout=zg(negvals_nh), method=methodzerzg, triangles=tr)
         gg(negvals_nh)=GRIDDATA(xp, zp, g, xout=xg(negvals_nh), yout=zg(negvals_nh), method=methodzerzg, triangles=tr)
         fhg(negvals_nh)=GRIDDATA(xp, zp, fh, xout=xg(negvals_nh), yout=zg(negvals_nh), method=methodzerzg, triangles=tr)
         tgg(negvals_nh)=GRIDDATA(xp, zp, tg, xout=xg(negvals_nh), yout=zg(negvals_nh), method=methodzerzg, triangles=tr)
         fwaterg(negvals_nh)=GRIDDATA(xp, zp, fwater, xout=xg(negvals_nh), yout=zg(negvals_nh), method=methodzerzg, triangles=tr)
         fohg(negvals_nh)=GRIDDATA(xp, zp, foh, xout=xg(negvals_nh), yout=zg(negvals_nh), method=methodzerzg, triangles=tr)
         
       endif
       
       ;zervals_nh = where(atan(zg/xg)*180.0/!pi le 3.0)
       zervals_nh = where(zg eq 0)
       if ((zervals_nh(0) ne -1)) then begin
         nhg(zervals_nh)=GRIDDATA(xpII, zpII, nhgII, xout=xg(zervals_nh), yout=zg(zervals_nh), method=methodzerzg, triangles=trII)
         epsig(zervals_nh)=GRIDDATA(xp, zp, epsi, xout=xg(zervals_nh), yout=zg(zervals_nh), method=methodzerzg, triangles=tr)
         scag(zervals_nh)=GRIDDATA(xp, zp, sca, xout=xg(zervals_nh), yout=zg(zervals_nh), method=methodzerzg, triangles=tr)
         absg(zervals_nh)=GRIDDATA(xp, zp, abss, xout=xg(zervals_nh), yout=zg(zervals_nh), method=methodzerzg, triangles=tr)
         gg(zervals_nh)=GRIDDATA(xp, zp, g, xout=xg(zervals_nh), yout=zg(zervals_nh), method=methodzerzg, triangles=tr)
         fhg(zervals_nh)=GRIDDATA(xp, zp, fh, xout=xg(zervals_nh), yout=zg(zervals_nh), method=methodzerzg, triangles=tr)
         tgg(zervals_nh)=GRIDDATA(xp, zp, tg, xout=xg(zervals_nh), yout=zg(zervals_nh), method=methodzerzg, triangles=tr)
         fwaterg(zervals_nh)=GRIDDATA(xp, zp, fwater, xout=xg(zervals_nh), yout=zg(zervals_nh), method=methodzerzg, triangles=tr)
         fohg(zervals_nh)=GRIDDATA(xp, zp, foh, xout=xg(zervals_nh), yout=zg(zervals_nh), method=methodzerzg, triangles=tr)
       endif


       methodlabelold = methodlabel 
       methodlabel = 'NaturalNeighbor'
       ; Fix the outer-most R values, had gotten weird interpolations.
       r_outermost = where(xg ge grid.xupper*0.95)
       ;print, xg(r_outermost)
       if ((r_outermost(0) ne -1)) then begin
         nhg(r_outermost)=GRIDDATA(xpII, zpII, nhgII, xout=xg(r_outermost), yout=zg(r_outermost), method=methodzerzg, triangles=trII)
         ;nhg(r_outermost) = 1e-12
         epsig(r_outermost)=GRIDDATA(xp, zp, epsi, xout=xg(r_outermost), yout=zg(r_outermost), method=methodzerzg, triangles=tr)
         scag(r_outermost)=GRIDDATA(xp, zp, sca, xout=xg(r_outermost), yout=zg(r_outermost), method=methodzerzg, triangles=tr)
         absg(r_outermost)=GRIDDATA(xp, zp, abss, xout=xg(r_outermost), yout=zg(r_outermost), method=methodzerzg, triangles=tr)
         gg(r_outermost)=GRIDDATA(xp, zp, g, xout=xg(r_outermost), yout=zg(r_outermost), method=methodzerzg, triangles=tr)
         fhg(r_outermost)=GRIDDATA(xp, zp, fh, xout=xg(r_outermost), yout=zg(r_outermost), method=methodzerzg, triangles=tr)
         tgg(r_outermost)=GRIDDATA(xp, zp, tg, xout=xg(r_outermost), yout=zg(r_outermost), method=methodzerzg, triangles=tr)
         fwaterg(r_outermost)=GRIDDATA(xp, zp, fwater, xout=xg(r_outermost), yout=zg(r_outermost), method=methodzerzg, triangles=tr)
         fohg(r_outermost)=GRIDDATA(xp, zp, foh, xout=xg(r_outermost), yout=zg(r_outermost), method=methodzerzg, triangles=tr)
       endif
       methodlabel = methodlabelold
       ;zervals_nh = where(zg eq 0)
       ;if ((zervals_nh(0) ne -1)) then begin
       ;  nhg(zervals_nh) = 0.0
       ;endif
       
       endif else begin
       
;      		for i=0,n_elements(xg)-1 do begin
;      		
;        		rsq=(xp-xg[i])^2.+(zp-zg[i])^2.
;            if(xg[i] lt grid.xgap) then begin
;              nhg[i]=0.
;            endif else begin
;              nhg[i]=nh[where(rsq EQ min(rsq))]
;            endelse
;      
;      		endfor
;      
       endelse

    ;=========== check for inside gap
      
    whingap=where(xg LT grid.xgap, numgap)
    if(numgap gt 0) then  begin
    
      nhg[whingap]=0.
      epsig[whingap]=0.
      scag[whingap]=0.
      absg[whingap]=0.
      gg[whingap]=0.
      fhg[whingap]=0.
      tgg[whingap]=1000.
      fwaterg[whingap]=0.
      fohg[whingap]=0.
      
    endif
      
    ;=========== finally insert valid grid points into node list
    
        xp=[xp,xg]
        zp=[zp,zg]
        
        nh=[nh,nhg] 
        epsi=[epsi,epsig]
        sca=[sca,scag]
        abss=[abss,absg]
        g=[g,gg]
        fh=[fh,fhg]
        tg=[tg,tgg]
        fwater=[fwater,fwaterg]
        foh=[foh,fohg]
        
        
        badvals_nh = where(finite(nh) eq 0.0)
        negvals_nh = where(nh lt 0.0)
        ;if badvals_nh(0) ne -1 or negvals_nh(0) ne -1 then begin
        ;  print,'w'
        ;endif
endif

methodzerzg = 'NaturalNeighbor'

;===========================================================
;  Cast a power-law radial grid, so cells have maximum size, helps plotting. 
;===========================================================

if(flags.castradial EQ 1) then begin

  a_castradial, xr, zr
  
    print, 'assigning densities to grid points...'
  
  nhr=0.*xr
 
 if(dotrianginterp EQ 1) then begin

       triangulate, xp, zp, tr

       zervec_xp = findgen(100)/101.0
       zervec_zp = fltarr(100)
       zervec_nh = dblarr(100)
       
       for jz = 0,99 do begin
         junk = min(sqrt((xp-zervec_xp(jz))^2+(zp-zervec_zp(jz))^2),minind)
         zervec_nh(jz) = nh(minind)
       endfor
       
       zervec_zp = fltarr(100)-0.05
       nhgII = [nh,zervec_nh]
       xpII =  [xp,zervec_xp]
       zpII =  [zp,zervec_zp]

       triangulate, xpII, zpII, trII
       
       nhr=GRIDDATA(xpII, zpII, nhgII, xout=xr, yout=zr, method=methodlabel, triangles=trII)
       epsir=GRIDDATA(xp, zp, epsi, xout=xr, yout=zr, method=methodlabel, triangles=tr)
       scar=GRIDDATA(xp, zp, sca, xout=xr, yout=zr, method=methodlabel, triangles=tr)
       absr=GRIDDATA(xp, zp, abss, xout=xr, yout=zr, method=methodlabel, triangles=tr)
       gr=GRIDDATA(xp, zp, g, xout=xr, yout=zr, method=methodlabel, triangles=tr)
       fhr=GRIDDATA(xp, zp, fh, xout=xr, yout=zr, method=methodlabel, triangles=tr)
       tgr=GRIDDATA(xp, zp, tg, xout=xr, yout=zr, method=methodlabel, triangles=tr)
       fwaterr=GRIDDATA(xp, zp, fwater, xout=xr, yout=zr, method=methodlabel, triangles=tr)
       fohr=GRIDDATA(xp, zp, foh, xout=xr, yout=zr, method=methodlabel, triangles=tr)

       badvals_nhr = where(finite(nhr) eq 0.0)
       negvals_nhr = where(nhr lt 0.0)

       if negvals_nhr(0) ne -1 then begin
         nhr(negvals_nhr)=GRIDDATA(xpII, zpII, nhgII, xout=xr(negvals_nhr), yout=zr(negvals_nhr), method=methodzerzg, triangles=trII)
         epsir(negvals_nhr)=GRIDDATA(xp, zp, epsi, xout=xr(negvals_nhr), yout=zr(negvals_nhr), method=methodzerzg, triangles=tr)
         scar(negvals_nhr)=GRIDDATA(xp, zp, sca, xout=xr(negvals_nhr), yout=zr(negvals_nhr), method=methodzerzg, triangles=tr)
         absr(negvals_nhr)=GRIDDATA(xp, zp, abss, xout=xr(negvals_nhr), yout=zr(negvals_nhr), method=methodzerzg, triangles=tr)
         gr(negvals_nhr)=GRIDDATA(xp, zp, g, xout=xr(negvals_nhr), yout=zr(negvals_nhr), method=methodzerzg, triangles=tr)
         fhr(negvals_nhr)=GRIDDATA(xp, zp, fh, xout=xr(negvals_nhr), yout=zr(negvals_nhr), method=methodzerzg, triangles=tr)
         tgr(negvals_nhr)=GRIDDATA(xp, zp, tg, xout=xr(negvals_nhr), yout=zr(negvals_nhr), method=methodzerzg, triangles=tr)
         fwaterr(negvals_nhr)=GRIDDATA(xp, zp, fwater, xout=xr(negvals_nhr), yout=zr(negvals_nhr), method=methodzerzg, triangles=tr)
         fohr(negvals_nhr)=GRIDDATA(xp, zp, foh, xout=xr(negvals_nhr), yout=zr(negvals_nhr), method=methodzerzg, triangles=tr)
       endif
       
       zervals_nhr = where(zr eq 0.0)
       if zervals_nhr(0) ne -1 then begin
         nhr(zervals_nhr)=GRIDDATA(xpII, zpII, nhgII, xout=xr(zervals_nhr), yout=zr(zervals_nhr), method=methodzerzg, triangles=trII)
         epsir(zervals_nhr)=GRIDDATA(xp, zp, epsi, xout=xr(zervals_nhr), yout=zr(zervals_nhr), method=methodzerzg, triangles=tr)
         scar(zervals_nhr)=GRIDDATA(xp, zp, sca, xout=xr(zervals_nhr), yout=zr(zervals_nhr), method=methodzerzg, triangles=tr)
         absr(zervals_nhr)=GRIDDATA(xp, zp, abss, xout=xr(zervals_nhr), yout=zr(zervals_nhr), method=methodzerzg, triangles=tr)
         gr(zervals_nhr)=GRIDDATA(xp, zp, g, xout=xr(zervals_nhr), yout=zr(zervals_nhr), method=methodzerzg, triangles=tr)
         fhr(zervals_nhr)=GRIDDATA(xp, zp, fh, xout=xr(zervals_nhr), yout=zr(zervals_nhr), method=methodzerzg, triangles=tr)
         tgr(zervals_nhr)=GRIDDATA(xp, zp, tg, xout=xr(zervals_nhr), yout=zr(zervals_nhr), method=methodzerzg, triangles=tr)
         fwaterr(zervals_nhr)=GRIDDATA(xp, zp, fwater, xout=xr(zervals_nhr), yout=zr(zervals_nhr), method=methodzerzg, triangles=tr)
         fohr(zervals_nhr)=GRIDDATA(xp, zp, foh, xout=xr(zervals_nhr), yout=zr(zervals_nhr), method=methodzerzg, triangles=tr)
       endif

       methodlabelold = methodzerzg
       methodzerzg = 'NaturalNeighbor'
       
       ; Fix the outer-most R values, had gotten weird interpolations.
       r_outermost = where(xr ge grid.xupper*0.95)

       if ((r_outermost(0) ne -1)) then begin       
         nhr(r_outermost)=GRIDDATA(xpII, zpII, nhgII, xout=xr(r_outermost), yout=zr(r_outermost), method=methodzerzg, triangles=trII)
         ;nhr(r_outermost)=1e-12
         epsir(r_outermost)=GRIDDATA(xp, zp, epsi, xout=xr(r_outermost), yout=zr(r_outermost), method=methodzerzg, triangles=tr)
         scar(r_outermost)=GRIDDATA(xp, zp, sca, xout=xr(r_outermost), yout=zr(r_outermost), method=methodzerzg, triangles=tr)
         absr(r_outermost)=GRIDDATA(xp, zp, abss, xout=xr(r_outermost), yout=zr(r_outermost), method=methodzerzg, triangles=tr)
         gr(r_outermost)=GRIDDATA(xp, zp, g, xout=xr(r_outermost), yout=zr(r_outermost), method=methodzerzg, triangles=tr)
         fhr(r_outermost)=GRIDDATA(xp, zp, fh, xout=xr(r_outermost), yout=zr(r_outermost), method=methodzerzg, triangles=tr)
         tgr(r_outermost)=GRIDDATA(xp, zp, tg, xout=xr(r_outermost), yout=zr(r_outermost), method=methodzerzg, triangles=tr)
         fwaterr(r_outermost)=GRIDDATA(xp, zp, fwater, xout=xr(r_outermost), yout=zr(r_outermost), method=methodzerzg, triangles=tr)
         fohr(r_outermost)=GRIDDATA(xp, zp, foh, xout=xr(r_outermost), yout=zr(r_outermost), method=methodzerzg, triangles=tr)
       endif
       
       methodzerzg = methodlabelold

 endif else begin

;    for i=0,n_elements(xr)-1 do begin
;    
;      rsq=(xp-xr[i])^2.+(zp-zr[i])^2.
;      if(xr[i] lt grid.xgap) then begin
;        nhr[i]=0.
;      endif else begin
;        nhr[i]=nh[where(rsq EQ min(rsq))]
;      endelse
;
;    endfor
    
 endelse

    ;=========== check for inside gap
      
    whingap=where(xr LT grid.xgap, numgap)
    if(numgap gt 0) then begin
    
      nhr[whingap]=0.0
      epsir[whingap]=0.
      scar[whingap]=0.
      absr[whingap]=0.
      gr[whingap]=0.
      fhr[whingap]=0.
      tgr[whingap]=1000.  ; prevent Doppler Nan and Infs
      fwaterr[whingap]=0.
      fohr[whingap]=0.
       
    endif
      
    ;=========== finally insert valid grid points into node list
    
        xp=[xp,xr]
        zp=[zp,zr]
        
        nh=[nh,nhr] 
        epsi=[epsi,epsir]
        sca=[sca,scar]
        abss=[abss,absr]
        g=[g,gr]
        fh=[fh,fhr]
        tg=[tg,tgr]
        fwater=[fwater,fwaterr]
        foh=[foh,fohr]
        whingapII=where(xp LT grid.xgap, numgap)
        if(numgap gt 0) then begin
          nh[whingapII]=0.0
          epsi[whingapII]=0.
          sca[whingapII]=0.
          abss[whingapII]=0.
          g[whingapII]=0.
          fh[whingapII]=0.
          tg[whingapII]=1000.  ; prevent Doppler Nan and Infs
          fwater[whingapII]=0.
          foh[whingapII]=0.
        endif
        
endif

;===========================================================
;  Cast corners if no background grids used at all
;===========================================================

;if(flags.castgrid + flags.castradial EQ 0) then begin
if(flags.castgrid EQ 0) then begin

    xp=[xp,0,0,grid.xupper,grid.xupper]
    zp=[zp,0,grid.zupper,0,grid.zupper]
    nh=[nh,0,0,0,0]
    epsi=[epsi,0,0,0,0]
    sca=[sca,0,0,0,0]
    abss=[abss,0,0,0,0]
    g=[g,0.6,0.6,0.6,0.6]
    fh=[fh,0,0,0,0]
    tg=[tg,1000,1000,1000,1000]
    fwater=[fwater,0,0,0,0]
    foh=[foh,0,0,0,0]

endif

;xp=xpl
;zp=zpl

;===========================================================
;  Check for point dupication.  Nudge.
;===========================================================

nudges=0L
seed = !NULL ; Initialize random seed for nudges in random +/- direction by 1.001x
print, nudges,' Testing for duplicated nodes...' 

for i=0L,LONG(n_elements(xp)-1) do begin

    ;whdup=where(xp EQ xp[i],numdup)
    whdup=where((abs(zp-zp[i]) lt 1e-15) and (abs(xp-xp[i]) lt 1e-15),numdup)     ; Check both X and Z.

     if(numdup GT 1) then begin
        ;print, 'dupicate found.  Nudging.'
        ;print, xp[i],zp[i]
        nudges++
        xp[i] = xp[i] * (1.0 + 0.001 * signum(RANDOMU(seed)-0.5))
        zp[i] = zp[i] * (1.0 + 0.001 * signum(RANDOMU(seed)-0.5))
     endif
    
endfor
 
   print, nudges,' points were duplicated, so nudged.' 
; stop
end
