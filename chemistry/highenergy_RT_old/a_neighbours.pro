;pro a_neighbours, trang,trangx,trangy,ntr,trneigh,edgetype,domminx,xupper,yupper
  pro a_neighbours

common GlobalVars

;====================================================
; Given Delaunay cell vertices.  turn these into quasi-unique edge values somehow, 
; log +log each vertex pari?  Method to so much important as spreading out these values 
; within the range of the numerical precision available.  It's a checkable chance method,
; but recall the approximate nearest neighbour routines...
;====================================================

	;trang=trang+!pi    ; just so that all elements >1

	tredid=dblarr(3,tr.ntr)
	tredid(0,*)=double(alog(tr.arr(0,*)+!pi))*double(alog(tr.arr(1,*)+!pi))  ; each edge gets unique value
	tredid(1,*)=double(alog(tr.arr(1,*)+!pi))*double(alog(tr.arr(2,*)+!pi))
	tredid(2,*)=double(alog(tr.arr(2,*)+!pi))*double(alog(tr.arr(0,*)+!pi))

	numtred=3.*tr.ntr  ;=3*ntr

;================================================================
; Identify edges along midplane y=0 etc for both vertices
;================================================================

	edgetypevalx=fltarr(3,tr.ntr)
	edgetypevalz=fltarr(3,tr.ntr)

	edgetypevalx(0,*)=tr.x(0,*)+tr.x(1,*)  ; combinations of coords that will characterize type of edge
	edgetypevalx(1,*)=tr.x(1,*)+tr.x(2,*)
	edgetypevalx(2,*)=tr.x(2,*)+tr.x(0,*)

	edgetypevalz(0,*)=tr.z(0,*)+tr.z(1,*)
	edgetypevalz(1,*)=tr.z(1,*)+tr.z(2,*)
	edgetypevalz(2,*)=tr.z(2,*)+tr.z(0,*)

	;==================================================
	;could change this so types take values 2,4,8,16 etc...
	; so dissimilar combinations are equally easily identified, c.f IDL style.
	;==================================================

	edgetype= fltarr(3,tr.ntr)   ; start with zeros
	edgetype=edgetype+2*(edgetypevalz eq 2.*grid.zupper)	; 2=top of domain
	edgetype=edgetype+3*(edgetypevalx eq 2.*grid.xupper)	; 3=Outer edge of domain
	edgetype=edgetype+4*(edgetypevalx eq 0.); 4=Inner edge of domain
	edgetype=edgetype+5*(edgetypevalz eq 0.)	; 5=midplane

tr.edgetype=edgetype

;================================================================
; turn into vectors, sort.
;================================================================

edgvec=reform(tredid,3.*tr.ntr)
celvec=LONG(findgen(LONG(3.*tr.ntr))/3.)     ; =(0,0,0,1,1,1,2,2,2, etc cell vector associated with tmpar
couvec=LONG(findgen(LONG(3.*tr.ntr)) mod 3)  ; =(0,1,2,0,1,2,0,1,2,vector indexing edges in systematic order

edgsort=LONG(sort(edgvec, /L64))   ; elemental sort order

edgvecs=edgvec[edgsort]    ; sorted into edge value order
celvecs=celvec[edgsort]   ; cells associated with above order
couvecs=couvec[edgsort]

;print, 'edgvec  :',edgvec
;print, 'celvec  :',celvec
;print, 'Now sorted...'
;print, 'edgvecs :',edgvecs
;print, 'celvecs :',celvecs


;================================================================
;
;================================================================

;tr.neigh=LONG(intarr(3,ntr))
tr.neigh(*,*)=-1L    ; start with -1.  

nedgvec=n_elements(edgvecs)
edgcount=LONG(intarr(tr.ntr))

print, 'Starting loop over edges'

for i=0L, LONG(nedgvec-2) do begin   ; loop over edges, ordered by their fictitious values.

	if(edgvecs(i) eq edgvecs(i+1)) then begin

		;print, ''
		;print, cellhere, cellnext, 'neighbors'
		;print, couvecs(cellhere),couvecs(cellnext)

		cellhere=LONG(celvecs(i))
		cellnext=LONG(celvecs(i+1))
		
		tr.neigh(couvecs(i),cellhere)=cellnext  ;this cell
		tr.neigh(couvecs(i+1),cellnext)=cellhere  ;complementary cell

	endif

endfor

end
