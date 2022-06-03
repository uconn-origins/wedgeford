pro a_castradial, xr,zr

common GlobalVars

print, 'Casting radial grid. nx,ny :',grid.nxr,grid.nzr
print, '... with node density r^-p, p=', flags.castradialpowerx

  ntot=grid.nxr*grid.nzr
  xr=dblarr(ntot)
  zr=dblarr(ntot)

  pts=0L


radius=sqrt(grid.xupper*grid.xupper+grid.zupper*grid.zupper)
rmax=1.414*radius
rmin=grid.xgap
rmaxminterm=(rmax/rmin)^(1.-flags.castradialpowerx)-1.

  for j=0,grid.nzr-1 do begin   ; altitude
    ;ang=j*0.5*!pi/(grid.nzr-1.)
    ang=j*(90.0-para.theta1)*!pi/180.0*1.25/(grid.nzr-1.)
    ;ang=
    for i=0,grid.nxr-1 do begin   ; radial
  
      wibble=(1.+0.001*(randomu(seed,1)-0.5))
      
  radi=rmin*(rmaxminterm*(1.*i/grid.nxr)+1.0)^(1./(1.-flags.castradialpowerx))
   xr[pts]=radi*cos(ang) 
   zr[pts]=radi*sin(ang) 
      ;
      ;xr(pts)=wibble*cos(ang)*1.414*radius*(0.155+0.845*i/(grid.nxg-1.))
      ;zr(pts)=wibble*sin(ang)*1.414*radius*(0.155+0.845*i/(grid.nxg-1.))
  
      pts=pts+1L
  
    endfor
  endfor


 whgood=where(xr LE grid.xupper AND zr LE grid.zupper,numgood)
  
  print, 'Of these ', numgood, ' fall in specified domain.'
  
  xr=xr[whgood]
  zr=zr[whgood]
 
  ;xr(0)=0.D
  ;zr(0)=0.D

;======== don't forget corners - DONE IN CASTINGMAIN

 ;  xr=[xr,0,0,grid.xupper,grid.xupper]
 ;  zr=[zr,0,grid.zupper,0,grid.zupper]

end