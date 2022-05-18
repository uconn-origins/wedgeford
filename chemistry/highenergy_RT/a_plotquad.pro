pro a_plotquad;, winny;,xplotmin,xplotscl,yplotscl

common GlobalVars

  a=findgen(17)*(!pi*2/16.)
  usersym, 1.4*cos(a), 1.4*sin(a), /fill
  
  lmaxd=max(alog10(cell.nh))
 lmind=2.
 
 lminfh=-5
 lmaxfh=0.1
 
  dlevs=10.^(lmind+(lmaxd-lmind)*findgen(41)/40.)
 elevs=10.^(-4+4.0*findgen(41)/40.)
 
 fhlevs=10.^(lminfh+(lmaxfh-lminfh)*findgen(41)/40.)
 
 ;===========================================================
;  contours maps of various quantities like epsi, fh
;===========================================================
 winny=2
 window, winny, xs=600,ys=900
  !p.multi=[0,2,3]
  !p.charthick=1
  loadct, 39
;
     triangulate, cell.comx,cell.comz,tt
  contour, cell.nh, cell.comx*grid.xphys, cell.comz*grid.xphys, triangulation=tt, /fill, levels=dlevs,$
    xtitle='r [AU]', ytitle='z [AU]', /iso, charsize=2, title='nH',$
    xrange=[0,grid.xupper*grid.xphys], yrange=[0,grid.zupper*grid.xphys], /xs, /ys
  oplot, vert.x*grid.xphys,vert.z*grid.xphys, psym=3, color=255                   ; overplot generator points
 
   contour, cell.epsi, cell.comx*grid.xphys, cell.comz*grid.xphys, triangulation=tt, /fill, levels=elevs,$
    xtitle='r [AU]', ytitle='z [AU]', /iso, charsize=2, title='e',$
    xrange=[0,grid.xupper*grid.xphys], yrange=[0,grid.zupper*grid.xphys], /xs, /ys
  oplot, vert.x*grid.xphys,vert.z*grid.xphys, psym=3, color=255                   ; overplot generator points
 
   contour, cell.fh, cell.comx*grid.xphys, cell.comz*grid.xphys, triangulation=tt, /fill, levels=fhlevs,$
    xtitle='r [AU]', ytitle='z [AU]', /iso, charsize=2, title='f(H)',$
    xrange=[0,grid.xupper*grid.xphys], yrange=[0,grid.zupper*grid.xphys], /xs, /ys
  oplot, vert.x*grid.xphys,vert.z*grid.xphys, psym=3, color=255                   ; overplot generator points
 
;   contour, cell.nh, cell.comx*grid.xphys, cell.comz*grid.xphys, triangulation=tt, /fill, levels=dlevs,$
;    xtitle='r [AU]', ytitle='z [AU]', /iso, charsize=2, title='nH',$
;    xrange=[0,grid.xupper*grid.xphys], yrange=[0,grid.zupper*grid.xphys], /xs, /ys
;  oplot, vert.x*grid.xphys,vert.z*grid.xphys, psym=3, color=255                   ; overplot generator points
; 
;   contour, cell.nh, cell.comx*grid.xphys, cell.comz*grid.xphys, triangulation=tt, /fill, levels=dlevs,$
;    xtitle='r [AU]', ytitle='z [AU]', /iso, charsize=2, title='nH',$
;    xrange=[0,grid.xupper*grid.xphys], yrange=[0,grid.zupper*grid.xphys], /xs, /ys
;  oplot, vert.x*grid.xphys,vert.z*grid.xphys, psym=3, color=255                   ; overplot generator points
; 
;   contour, cell.nh, cell.comx*grid.xphys, cell.comz*grid.xphys, triangulation=tt, /fill, levels=dlevs,$
;    xtitle='r [AU]', ytitle='z [AU]', /iso, charsize=2, title='nH',$
;    xrange=[0,grid.xupper*grid.xphys], yrange=[0,grid.zupper*grid.xphys], /xs, /ys
;  oplot, vert.x*grid.xphys,vert.z*grid.xphys, psym=3, color=255                   ; overplot generator points
 
   
;===========================================================
;  10% innermost sub-region
;===========================================================
 winny=1
 window, winny, xs=220./(grid.zupper/grid.xupper),ys=960
  !p.multi=[0,1,4]
  !p.charthick=1
  loadct, 39
;
     triangulate, cell.comx,cell.comz,tt
  contour, cell.nh, cell.comx*grid.xphys, cell.comz*grid.xphys, triangulation=tt, /fill, levels=dlevs,$
    xtitle='r [AU]', ytitle='z [AU]', /iso, charsize=2, title='innermost 10%',$
    xrange=[0,0.1*grid.xupper*grid.xphys], yrange=[0,0.1*grid.zupper*grid.xphys], /xs, /ys
 
 oplot, vert.x*grid.xphys,vert.z*grid.xphys, psym=3, color=255                   ; overplot generator points
 
;===========================================================
;  1% innermost sub-region
;===========================================================
; winny=2
; window, winny, xs=220./(grid.zupper/grid.xupper),ys=220
;  !p.multi=[0,1,1]
;  !p.charthick=1.
;  loadct, 39
;
  triangulate, cell.comx,cell.comz,tt
  contour, cell.nh, cell.comx*grid.xphys, cell.comz*grid.xphys, triangulation=tt, /fill, levels=dlevs,$
    xtitle='r [AU]', ytitle='z [AU]', /iso, charsize=2, title='innermost 1%',$
    xrange=[0,0.01*grid.xupper*grid.xphys], yrange=[0,0.01*grid.zupper*grid.xphys], /xs, /ys
 
 oplot, vert.x*grid.xphys,vert.z*grid.xphys, psym=3, color=255                   ; overplot generator points

;===========================================================
;  0.1% innermost sub-region
;===========================================================
; winny=2
; window, winny, xs=220./(grid.zupper/grid.xupper),ys=220
;  !p.multi=[0,1,1]
;  !p.charthick=1.
;  loadct, 39
;
  triangulate, cell.comx,cell.comz,tt
  contour, cell.nh, cell.comx*grid.xphys, cell.comz*grid.xphys, triangulation=tt, /fill, levels=dlevs,$
    xtitle='r [AU]', ytitle='z [AU]', /iso, charsize=2, title='innermost 0.1%',$
     xrange=[0,0.001*grid.xupper*grid.xphys], yrange=[0,0.001*grid.zupper*grid.xphys], /xs, /ys
 
; for i=0L,tr.ntr-1L do begin
;    t=[tr.arr[*,i],tr.arr[0,i]]
;    plots, vert.x[t]*grid.xphys,vert.z[t]*grid.xphys, color=155
;    ;plots, vert.x[t]*grid.xphys,vert.z[t]*grid.xphys, color=0, thick=2
;    endfor
; 
 
 oplot, vert.x*grid.xphys,vert.z*grid.xphys, psym=3, color=255                   ; overplot generator points

;===========================================================
;  0.01% innermost sub-region
;===========================================================
; winny=2
; window, winny, xs=220./(grid.zupper/grid.xupper),ys=220
;  !p.multi=[0,1,1]
;  !p.charthick=1.
;  loadct, 39
;
  triangulate, cell.comx,cell.comz,tt
  contour, cell.nh, cell.comx*grid.xphys, cell.comz*grid.xphys, triangulation=tt, /fill, levels=dlevs,$
    xtitle='r [AU]', ytitle='z [AU]', /iso, charsize=2, title='innermost 0.01%',$
     xrange=[0,0.0001*grid.xupper*grid.xphys], yrange=[0,0.0001*grid.zupper*grid.xphys], /xs, /ys
 
; for i=0L,tr.ntr-1L do begin
;    t=[tr.arr[*,i],tr.arr[0,i]]
;    plots, vert.x[t]*grid.xphys,vert.z[t]*grid.xphys, color=155
;    ;plots, vert.x[t]*grid.xphys,vert.z[t]*grid.xphys, color=0, thick=2
;    endfor
 
 oplot, vert.x*grid.xphys,vert.z*grid.xphys, psym=3, color=255                   ; overplot generator points


;===========================================================
;  100% region
;===========================================================

winny=0
	window, winny, xs=720./(grid.zupper/grid.xupper),ys=720
	!p.multi=[0,1,1]
  !p.charthick=2.
	
;===========================================================
;  Plot density, generator points, triangulation.  One quadrant
;===========================================================

	loadct, 39
;  loadct, 1
 
  
 ; triangulate, cell.comx,cell.comz,tt
  contour, cell.nh, cell.comx*grid.xphys, cell.comz*grid.xphys, triangulation=tt, /fill, levels=dlevs,$
    xtitle='r [AU]', ytitle='z [AU]', /iso, charsize=2.4, xrange=[0,grid.xupper*grid.xphys], yrange=[0,grid.zupper*grid.xphys], /xs, /ys

;xyouts,cell.comx*grid.xphys, cell.comz*grid.xphys,strtrim(string(tr.edgetype[]))
	
;===========================================================
; k
;===========================================================

	loadct, 0, /silent

;*****
	for i=0L,tr.ntr-1L do begin
		t=[tr.arr[*,i],tr.arr[0,i]]
		plots, vert.x[t]*grid.xphys,vert.z[t]*grid.xphys, color=155
    ;plots, vert.x[t]*grid.xphys,vert.z[t]*grid.xphys, color=0, thick=2
  	endfor

	;oplot, vert.x*grid.xphys,vert.z*grid.xphys, psym=8, color=255                   ; overplot generator points
  oplot, vert.x*grid.xphys,vert.z*grid.xphys, psym=3, color=255                   ; overplot generator points
  ;oplot, orig.xp*grid.xphys,orig.zp*grid.xphys, psym=3, color=255                   ; overplot generator points

	;oplot, cell.comx*grid.xphys, cell.comz*grid.xphys, psym=3

;===========================================================
;  envelope into which we fire our photons
;===========================================================

oplot, [0,grid.xupper]*grid.xphys,[0,grid.xupper]*grid.xphys/tan(const.degtorad*para.theta1), color=254, linestyle=2, thick=3
oplot, [0,grid.xupper]*grid.xphys,[0,grid.xupper]*grid.xphys/tan(const.degtorad*para.theta2), color=30, linestyle=2, thick=3

;===========================================================
;  Plot star
;===========================================================

;  ang=!pi*0.5*findgen(51)/50.
;  oplot, 0.005*sin(ang)*grid.xphys,0.005*cos(ang)*grid.xphys, thick=4, color=200

;===========================================================
; "fuck the colourbar"
;===========================================================

;	colorbar,   maxrange=max(ldenlevs), minrange=min(ldenlevs), position=[0.92,0.1,0.96,0.9], /VERTICAL,$
	;		title='Log density', charsize=1.2, /RIGHT

;stop
  loadct, 39, /silent

end
