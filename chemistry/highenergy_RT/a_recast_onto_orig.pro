pro a_recast_onto_orig

common GlobalVars

;restore, './Output/testrun2_testdisk2_Jan_14_2010_blocks.sav'
  
  print, 'Recasting results onto the original set of disk coords ...'
  
  triangulate, cell.comx, cell.comz,tr
  orig.nph=GRIDDATA(cell.comx, cell.comz, ph.nph, xout=orig.xp, yout=orig.zp, method='linear', triangles=tr)
  print, 'here1'
  orig.fx=GRIDDATA(cell.comx, cell.comz, ph.nph, xout=orig.xp,yout=orig.zp, method='linear', triangles=tr)
  print, 'here1.5'
  orig.fz=GRIDDATA(cell.comx, cell.comz, ph.nph, xout=orig.xp,yout=orig.zp, method='linear', triangles=tr)

   print, 'here2'

;================================
; test these against gridded results
;================================

if(plotFlags.nodevice EQ 0) then begin
  if(plotFlags.origgrid) then begin
      
      window, 5, ys=400, xs=600
      !p.multi=[0,1,1]
      
        nl=101
      
        nhlevels=10.^(8.5*findgen(nl)/(nl-1.))
        nhlevels2=10.^(8.5*findgen(10)/(9.))
       
        jlevels=10.^(-4.*(1.-findgen(nl)/(nl-1.)))
        nphlevels=10.^(-1.+(8.*findgen(nl)/(nl-1.)))  
        alevels=1.*(findgen(nl)/(nl-1.))
        flevels=10.^(4.*findgen(nl)/(nl-1.))
      
      contour, orig.nph, orig.xp*grid.xphys, orig.zp*grid.xphys,/irregular,levels=nphlevels, /fill, $
          title='n!lph!n (recast onto orig grid)', /iso, xtitle='r [AU]', ytitle='z [AU]', xrange=[0, grid.xphys*grid.xupper], yrange=[0, grid.xphys*grid.zupper], /xs, /ys
   ;   contour, ph.nph,cell.comx*grid.xphys,cell.comz*grid.xphys,/irregular,levels=nphlevels, /fill, $
   ;       title='n!lH!n (orig grid)', /iso, xtitle='r [AU]', ytitle='z [AU]', xrange=[0, grid.xphys*grid.xupper], yrange=[0, grid.xphys*grid.zupper], /xs, /ys
   
   endif
 endif

end
