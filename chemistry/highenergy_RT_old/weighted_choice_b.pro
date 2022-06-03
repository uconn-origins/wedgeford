function weighted_choice_b,randomnum

  numelm = 500
  x = findgen(numelm)/float(numelm-1) * 3.14159
  weights = 1.0+(cos(x))^2.0
   
  totalsav = make_array(numelm,/float)

  running_total = 0.0
  

  for wc=0,numelm-1 do begin
    running_total = running_total + weights(wc)
    totalsav(wc) = running_total
  endfor
  

  rnd = randomnum * running_total
  
  minv = min(abs(totalsav-rnd),index)
  outval = float(index)/float(numelm) * 3.14159
  
  return, outval

end