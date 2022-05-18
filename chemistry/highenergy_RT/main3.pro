pro MAIN3;,lyphdgrid,waterfracgrid,dgrid,xax,zax

common GlobalVars

;==================================
;  


;==================================
; Plot and analyze results.  Plot unstructured, regrid. Zoom in
;==================================

;nl=101
;
;	nhlevels=10.^(8.5*findgen(nl)/(nl-1.))
;	nhlevels2=10.^(8.5*findgen(10)/(9.))
; 
; 	jlevels=10.^(-4.*(1.-findgen(nl)/(nl-1.)))
;	phdlevels=10.^(-1.+(8.*findgen(nl)/(nl-1.)))	
;  alevels=1.*(findgen(nl)/(nl-1.))
;  flevels=10.^(4.*findgen(nl)/(nl-1.))
; levs={numl:101,$
;        numlfew:11,$
;        nh:10.^(8.5*findgen(101)/100.),$
;        nhfew:10.^(8.5*findgen(11)/10.),$
;        j:10.^(8.5*findgen(101)/100.),$
;        nph:10.^(8.5*findgen(101)/100.),$
;        a:10.^(8.5*findgen(101)/100.),$
;        f:10.^(8.5*findgen(101)/100.)}

;============= 
; gridify
;=============

  nxax=128
  nzax=128
    
print, 'Making regular gridded version of results... '

p_gridify, nhgrid, nphgrid, fxgrid, fzgrid, xax, zax, nxax, nzax

 out={$
  nhgrid:nhgrid,$
  nphgrid:nphgrid,$
  fxgrid:fxgrid,$
  fzgrid:fzgrid,$
  fgrid:nhgrid,$
  agrid:nhgrid,$
  nxax:nxax,$
  nzax:nzax,$
  xax:fltarr(nxax),$
  zax:fltarr(nzax)}
 
print, 'out'
help, out.nhgrid,out.nphgrid
print, 'min(out.nphgrid),max(out.nphgrid) :',min(out.nphgrid),max(out.nphgrid)
print, 'min(out.nhgrid),max(out.nhgrid)   :',min(out.nhgrid),max(out.nhgrid)


;==========================
; plot results in unstructured form?
;==========================

if(plotFlags.nodevice EQ 0) then begin
	
  if(plotFlags.irr eq 1 AND plotFlags.global EQ 1) then begin

	window, 2, xsize=1000, ysize=400
	wset, 2
	
	!p.multi=[0,2,1]
	!p.charthick=1.
	!p.charsize=1.4

	;===========================================
	; Density, unsturctured
	;===========================================

	loadct, 1, /silent
	contour, cell.nh,cell.comx*grid.xphys,cell.comz*grid.xphys,/irregular,levels=levs.nh, /fill, $
		title='n!lH!n (unstructured)', /iso, xtitle='r', ytitle='z'
	
	contour, cell.nh,cell.comx*grid.xphys,cell.comz*grid.xphys,/irregular, /overplot, levels=levs.nhfew

     colorbar,   maxrange=max(alog10(levs.nh)), minrange=min(alog10(levs.nh)), position=[0.12,0.9,0.46,0.95], $
      title='Log n!lH!n', charsize=1.2, /RIGHT

	plots, [0],[0],psym=2, thick=6

	;===========================================
	; photon number density unstructured
	;===========================================

	loadct, 39, /silent
  contour,ph.nph,cell.comx*grid.xphys,cell.comz*grid.xphys,/irregular, /fill,  /iso, $
		levels=levs.nph,title='n!lph!n (unstructured) ph cm!u-3!n', xtitle='r', ytitle='z'

     colorbar,   maxrange=max(alog10(levs.nph)), minrange=min(alog10(levs.nph)), position=[0.62,0.9,0.96,0.95], $
      title='Log n!lph!n', charsize=1.2, /RIGHT

  endif
endif

;===========================================
; Plot regularly gridded version
;===========================================

       out.fgrid=sqrt(out.fxgrid*out.fxgrid+out.fzgrid*out.fzgrid)   ; magnitude of flux
       out.agrid=out.fgrid/(out.nphgrid*3.E10)

       finiteagrid=finite(out.agrid)
       out.agrid(where(finiteagrid ne 1))=0.01
       
;       amask=0.*acgrid+1.
;       
;        for j=0,nzax-1 do begin
;          for i=0,nxax-1 do begin
;          
;            if(1.*j/i GE 3.) then amask[i,j]=0.01
;           
;          endfor
;        endfor
; 
;      agrid=agrid*amask

        fac=4
        print, 'fac =',fac
        
        vecfield=fltarr(2,nxax/fac,nzax/fac)
        vecfield(0,*,*)=rebin(out.fxgrid/out.fgrid,nxax/fac,nzax/fac)
        vecfield(1,*,*)=rebin(out.fzgrid/out.fgrid,nxax/fac,nzax/fac)
        rebx=grid.xupper*findgen(nxax/fac)/(nxax/fac-1.)
        rebz=grid.zupper*findgen(nzax/fac)/(nzax/fac-1.)
 
        vecx=fltarr(nxax/fac,nzax/fac)
        vecz=vecx
        vecx(*,*)=vecfield(0,*,*)
        vecz(*,*)=vecfield(1,*,*)
 
        
        ;if(plotFlags.reg eq 1) then begin
        
 if(plotFlags.nodevice EQ 0) then begin
 
  if(plotFlags.gridded eq 1 AND plotFlags.global EQ 1) then begin
 
         print, 'Plotting regular gridded version of results... '
        
        	window, 3, xsize=1000, ysize=800
        	wset, 3
        	
        	!p.multi=[0,2,2]
        	!p.charsize=1.4
        	!p.charthick=1.
        	
          ;============= Density on grid
  
        	loadct, 1, /silent
        
        	contour,out.nhgrid,xax*grid.xphys,zax*grid.xphys, /fill, /iso, $
        		title='n!lH!n (gridded)', xtitle='r', ytitle='z', levels=levs.nh
         
          contour,out.nhgrid,xax*grid.xphys,zax*grid.xphys, levels=levs.nhfew, /overplot

    colorbar,   maxrange=max(alog10(levs.nh)), minrange=min(alog10(levs.nh)), position=[0.12,0.9,0.46,0.95], $
    title='Log n!lH!n', charsize=1.2, /RIGHT
        
        ;============= Photon density on grid
   
        	loadct, 39, /silent
        	contour, out.nphgrid,xax*grid.xphys,zax*grid.xphys,levels=levs.nph, /fill, $
        		title='n!lph!n (gridded)', /iso, xtitle='r', ytitle='z'
 
          contour, out.nphgrid,xax*grid.xphys,zax*grid.xphys,levels=levs.nphfew, /overplot, $
            title='n!lph!n (gridded)', /iso, xtitle='r', ytitle='z', C_LABELS=replicate(1,10), C_CHARSIZE=1.2

    colorbar,   maxrange=max(alog10(levs.nph)), minrange=min(alog10(levs.nph)), position=[0.62,0.9,0.96,0.95], $
    title='Log n!lph!n', charsize=1.2, /RIGHT
         	
      	;====== anisotropy of flux ======
      ;
      ;  window, 4, xsize=1000, ysize=400
      ; 
      ;  !p.multi=[0,2,1]
      ;  !p.charthick=1.
      ;  !p.charsize=1.4

      contour,out.fgrid,xax*grid.xphys,zax*grid.xphys, /fill, /iso, $
      title='flux (gridded)', xtitle='r', ytitle='z', levels=levs.f
 
     colorbar,   maxrange=max(alog10(levs.f)), minrange=min(alog10(levs.f)), position=[0.12,0.4,0.5,0.45], $
      title='log10(Flux [ph cm-2])', charsize=1.2, /RIGHT

      contour,out.agrid,xax*grid.xphys,zax*grid.xphys, /fill, /iso, $
      title='anisotropy (gridded)', xtitle='r', ytitle='z', levels=levs.a
               
        velovect, vecx, vecz, rebx*grid.xphys,rebz*grid.xphys, /iso, /overplot, color=0, thick=1.4

      colorbar,   maxrange=1, minrange=0, position=[0.62,0.4,0.96,0.45], $
      title='anisotropy', charsize=1.2, /RIGHT

        ;endif
         
          ;=====================================
          ;  vertical slice: reality check
          ;=====================================
        
        ;window, 4
        ;!p.multi=[0,1,2]
        radsel=grid.xphys*grid.xupper*0.5
        whpts=where(cell.comx*grid.xphys GE radsel AND cell.comx*grid.xphys LE radsel*1.1, numwh)
        
        if(numwh GT 0) then begin
        
         window, 4, xs=400,ys=800, title='Reality check'
         !p.multi=[0,1,2]
          ;print, 'Vertical selection of ',numwh,' points.',ph.nph[whpts]
        
         ;========== photon density
          plot_IO, cell.comz[whpts]*grid.xphys, ph.nph[whpts], yrange=[1E-4,10]*max(ph.nph[whpts]), /ys, psym=1, ytitle='n_ph',$
           xtitle='Z [AU]', color=255, title='vertical sample at '+string(radsel)+'AU'
        
          rsq=((cell.comz[whpts]*grid.xphys)*(cell.comz[whpts]*grid.xphys)+(cell.comx[whpts]*grid.xphys)*(cell.comx[whpts]*grid.xphys))*2.25D26
          nphoton=para.starlumin/(4.*!pi*rsq*3.0D10)
        
         oplot, cell.comz[whpts]*grid.xphys, nphoton, psym=3, color=254
         
         ;========== flux
         
          fphoton=para.starlumin/(4.*!pi*rsq)
         
         plot_IO, cell.comz[whpts]*grid.xphys, sqrt(ph.fx[whpts]*ph.fx[whpts]+ph.fz[whpts]*ph.fz[whpts]), yrange=[1E-4,10]*max(sqrt(ph.fx[whpts]*ph.fx[whpts]+ph.fz[whpts]*ph.fz[whpts])), /ys, psym=1, ytitle='FLux', xtitle='Z [AU]', color=255
         oplot, cell.comz[whpts]*grid.xphys, fphoton, psym=3, color=254
        
        endif
        
          print, '---- unattenuated values (i.e. inverse square), based on stellar props ----'
          print, 'Star luminosity 1E44 phot/s yields:'
          print, 'nph cm-3=',1D44/(4.*!pi*2.25E26*3.0E10),' at distance of 1AU'
          print, 'nph cm-3=',1D44/(4.*!pi*10000*2.25D26*3.0D10),' at distance of 100AU'
          print, ' '
          print, 'Your star of luminosity ',para.starlumin,' phot/s yields:'
          print, 'nph cm-3=',para.starlumin/(4.*!pi*2.25E26*3.0E10),' at distance of 1AU'
          print, 'nph cm-3=',para.starlumin/(4.*!pi*10000*2.25D26*3.0D10),' at distance of 100AU'
          print, ' '
          print, 'Flux ph/cm2/s =',para.starlumin/(4.*!pi*2.25E26),' at distance of 1AU'
          print, 'Flux ph/cm2/s =',para.starlumin/(4.*!pi*10000*2.25D26),' at distance of 100AU'

 endif
endif
	;===================================================================
	;  
	;===================================================================




;stop
end
