pro a_disk_create, xp, zp, nhp, ntr, comx, comz, nh, epsi, sca, abss, g, fh, tg, fwater, foh, trang, trangx, trangz, trar, trvol
 
; xp,zp and nhp are the values at the vertices.  
; the vlues at cell centers (comx) are without p suffix.

;===========================================================
;  Load disk, cast background grids.
;  Returns values at nodes.
;===========================================================

a_casting_main, xp, zp, nhp, epsip, scap, absp, gp, fhp, tgp, fwaterp, fohp

print, 'Nodes finalized.  Triangulating ...'

;===========================================================
;  Triangulate.  Make array of vertex coords.  Density at vertices is known.
;===========================================================

	triangulate,xp,zp,trang,bounder,connectivity=tlist  ;trang is array of indices, ie LONG

	ntr=n_elements(trang)/3
	trangx=xp(trang)    ; trangx lists the x coord of vertices j connected to it.
	trangz=zp(trang)

print, ntr, ' cells generated in domain.'

;===========================================================
;  Areas, volumes, densities, linear scale, columns of cells.
;===========================================================

	print, 'Calculating areas, CoM and density of cells'

	a_areas_vols, trang, ntr, xp, zp, trlin, trar, trvol, comx, comz

;===========================================================
;  Interpolate disk data onto comx,comz
;  These are CELL CENTERS.
;===========================================================

nh=fltarr(ntr)
epsi=fltarr(ntr)
abss=fltarr(ntr)
sca=fltarr(ntr)
g=fltarr(ntr)
fh=fltarr(ntr)
tg=fltarr(ntr)
fwater=fltarr(ntr)
foh=fltarr(ntr)

    for i=0L,LONG(ntr)-1L do begin
    
        nh[i]=total(nhp[trang[*,i]])/3.
        epsi[i]=total(epsip[trang[*,i]])/3.
        abss[i]=total(absp[trang[*,i]])/3.
        sca[i]=total(scap[trang[*,i]])/3.
        g[i]=total(gp[trang[*,i]])/3.
        fh[i]=total(fhp[trang[*,i]])/3.
        tg[i]=total(tgp[trang[*,i]])/3.
        fwater[i]=total(fwaterp[trang[*,i]])/3.
        foh[i]=total(fohp[trang[*,i]])/3.
    
    endfor

end
