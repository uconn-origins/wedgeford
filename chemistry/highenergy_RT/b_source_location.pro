pro B_source_location, comx,comy, fixti

method=3

; =0 cell nearest origin
; =1 specify cell index
; =2 random cell
; =3 cell nearest specified point
;

if(method eq 2) then begin
	fixti=LONG(ntr*randomu(seed))  ; initial cell
endif

if(method eq 1) then begin
	fixti=510
endif

if(metod eq 0) then begin
	gz=sqrt((comx*comx)+(comz*comz))
	wgz=where(gz eq min(gz))
	fixti=wgz(0)
endif


if(method eq 3) then begin

	;  given chosen position, find nearest cell com and use that instead.

	chpos=[19.5,7.5]
	chpos=[0.,0.]

	diff=(comx-chpos(0))*(comx-chpos(0))+(comz-chpos(1))*(comz-chpos(1))
	fixti=where(diff eq min(diff))

endif


end
