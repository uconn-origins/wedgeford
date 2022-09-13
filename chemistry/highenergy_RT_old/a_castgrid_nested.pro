pro a_castgrid_nested, xg, zg
 
 common GlobalVars 
;common domainblock2
;common trangblock2
;common flagblock2

 	print, 'Casting regular grid. nx,ny :',grid.nxg,grid.nzg

	if(grid.nxg lt 2) then grid.nxg=2  ; ensures that four corners will be done
	if(grid.nzg lt 2) then grid.nzg=2

	ntot=grid.nxg*grid.nzg
	xg=dblarr(ntot)
	zg=dblarr(ntot)

	pts=0L

	xsample=double(grid.xupper)*dindgen(grid.nxg)/(DOUBLE(grid.nxg)-1.D)   ; over entire domain
	zsample=double(grid.zupper)*dindgen(grid.nzg)/(DOUBLE(grid.nzg)-1.D)

	for j=0,grid.nzg-1 do begin
	for i=0,grid.nxg-1 do begin

		xg(pts)=xsample(i)		;+0.01D*randomu(seed,1)
		zg(pts)=zsample(j)		;+0.01D*randomu(seed,1)

		pts=pts+1L

	endfor
	endfor

	xg(0)=0.D
	zg(0)=0.D

;===================================================
;  Cascading resolution in lower left corner.
; Half grid points in x,z, placed at midpoints
;===================================================

if(flags.castnested gt 0) then begin

  clevels=flags.castnested
  
  print, 'Adding a nested grid with ',clevels,' levels...'

for cl=1,clevels do begin

		xgc=dblarr(ntot)
		zgc=dblarr(ntot)
		pts=0L
		skipped=0L
		
		xsample=(0.5D^DOUBLE(cl))*double(grid.xupper)*dindgen(grid.nxg)/(double(grid.nxg)-1.D)   ; over entire domain
		zsample=(0.5D^DOUBLE(cl))*double(grid.zupper)*dindgen(grid.nzg)/(double(grid.nzg)-1.D)

	for j=0,grid.nzg-1 do begin
	for i=0,grid.nxg-1 do begin

		;if both even then skip as this point exists from previous level

		if(LONG(0.5*(i)) eq 0.5*(i) and LONG(0.5*(j)) eq 0.5*(j)) then begin

			skipped++
	
		endif else begin

			xgc(pts)=xsample(i)*(1.+0.01*randomu(seed,1))
			zgc(pts)=zsample(j)*(1.+0.01*randomu(seed,1))

			pts=pts+1L

		endelse

	endfor
	endfor

	xgc=[xgc(0:pts-1)]
	zgc=[zgc(0:pts-1)]

	;oplot, xgc,zgc, psym=3+cl

	xg=[xg,xgc]
	zg=[zg,zgc]

endfor

endif

;print, 'PRECISION OF GRID NODES :'
;help, xg,zg

;stop

end
