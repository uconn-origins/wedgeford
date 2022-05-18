
pro a_areas_vols,trang, ntr, x, y, trlin, trar, trvol, comx, comz

; Heron's formula A=sqrt(s(s-a)(s-b)(s-c))  where s is semi-perimeter 
; and a, b, and c are edge lengths.

trar=dblarr(ntr)
comx=dblarr(ntr)
comz=dblarr(ntr)
trlin=dblarr(ntr)

for i=0L,LONG(ntr-1) do begin

	x1=double(x(trang(0,i)))
	x2=double(x(trang(1,i)))
	x3=double(x(trang(2,i)))
	y1=double(y(trang(0,i)))
	y2=double(y(trang(1,i)))
	y3=double(y(trang(2,i)))

	al=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2))
	bl=sqrt((x2-x3)*(x2-x3)+(y2-y3)*(y2-y3))
	cl=sqrt((x3-x1)*(x3-x1)+(y3-y1)*(y3-y1))

	s=0.5D*(al+bl+cl)

	trar(i)=sqrt(s*(s-al)*(s-bl)*(s-cl))
	trlin(i)=sqrt((4.0D/!dpi)*trar(i))

	comx(i)=(x1+x2+x3)/3.0D
	comz(i)=(y1+y2+y3)/3.0D

endfor

trvol=trar*2.D*!dpi*comx  ;  No comz, rings are about x=0 axis.

print, 'Min cell linear scale :',sqrt(min(trar))
;print, 'Mean cell linear scale :',sqrt(mean(trar))
;stop

end
