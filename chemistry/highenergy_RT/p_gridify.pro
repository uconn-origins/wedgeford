pro p_gridify, nhgrid, nphgrid, fxgrid, fzgrid, xax, zax, nxax, nzax

common GlobalVars

;========================= 
;       gridify
;========================= 

	triangulate, cell.comx,cell.comz,triangul

  nphgrid=trigrid(cell.comx,cell.comz,ph.nph,triangul,nx=nxax,ny=nzax)
	nhgrid=trigrid(cell.comx,cell.comz,cell.nh,triangul,nx=nxax,ny=nzax)
  fxgrid=trigrid(cell.comx,cell.comz,ph.fx,triangul,nx=nxax,ny=nzax)
  fzgrid=trigrid(cell.comx,cell.comz,ph.fz,triangul,nx=nxax,ny=nzax)

	;xax=rmin+(rmax-rmin)*findgen(nxax)/(nxax-1.)
  xax=grid.xupper*findgen(nxax)/(nxax-1.)
	zax=grid.zupper*findgen(nzax)/(nzax-1.)

 ; help, nhgrid,nphgrid
 ; print, max(nhgrid),max(nphgrid)
 ; print, min(nhgrid),min(nphgrid)


end
