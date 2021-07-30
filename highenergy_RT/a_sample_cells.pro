pro A_sample_cells, comx, comz, dcell, sample, sampleel, nsample

;========================================
; find sample points by location
;========================================

	;sample=where(comx gt 0.8 and comx lt 0.81 and comz gt 0.1 and comz lt 0.2, nsample)

;========================================
; find sample points by density range 
;========================================

	lowerden=1.E5
	upperden=1.E9      ; replace this with blackden, or add that info?

	;sample=where(dcell ge lowerden and dcell le upperden and comx gt 0.2$
	;		 and comx gt 0.78 and comx lt 0.84, nsample)
	sample=where(dcell ge lowerden and dcell le upperden and comx gt 0.2, nsample)
	;sample=where(dcell ge lowerden and dcell le upperden and comx gt 0.8 and comx lt 0.21, nsample)

sample=[0,sample]
nsample++

;========================================
;
;========================================

	sampleel=0.*comx
	
	list=findgen(nsample)+1
	sampleel(sample)=list   ; for every cell, 0 not a sample cell, 
				; otherwise value = position+1 in samplehist array 

	print, 'Criterion -> sample of size :',nsample

	;oplot, comx(sample), comz(sample), psym=3, color=0	

end
