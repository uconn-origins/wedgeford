pro a_readparas, parafile;, parasf, parasl, parasc

common GlobalVars

;parasf=dblarr(30)
;parasc=strarr(30)
;parasl=lonarr(30)

un=3

prjdir='.'
savfile='saveas'
diskfile='diskread'

starlumin=0.0

    rmin=0.
    rmax=0.
    scrap='*'
    np=0
    nxg=0
    nzg=0
    castgrid=0
    castnested=0
    xfocus=0.0
    zfocus=0.0
    xupper=0.
    zupper=0.
    xphys=0.
  
    totnumph=1000000L
    theta1=0.
    theta2=0.
    wavelength=0.
    maxscats=0L
    lymanswitch=0
    randomizeangle=0
    usedustabs=1
    usedustsca=1
    xrayflag=0
    
    constabs=0.0E0
    constsca=0.0E0
    constepsilon=0.0000
    constg=0.0
   
   usemolabs=0.0
   
   plotrun=0
   plotirr=0
   plotreg=0
   
   
 ;   loadlyspecflag=0
 ;   loadlyspecfile=''
  ;  xmin=-1.
  ;  xmax=1.
  ;  wavmin=0.
  ;  wavmax=0.
  ;  numwav=0
  ;  restartitdomain=0L
  ;  restartitsource=0L
  ;  restartfilly=''
  ;  Mstar=0.
  ;  oabund=0.

;============================================
;
;============================================
;print, !PATH

print, 'Opening parameter file:', parafile

openr, un,parafile

readf, un, scrap   ;'12345678901234567890'
readf, un, savfile, FORMAT='(20X,A50)'
readf, un, scrap ;'-------------------- '
readf, un, diskfile, FORMAT='(20X,A50)'
readf, un, wavelength, FORMAT='(20X,F)'
readf, un, starlumin, FORMAT='(20X,D)'
readf, un, isrfg, FORMAT='(20X,D)'
readf, un, scrap ;'-------------------- Parameters for additional nested background grid and Poisson '
readf, un, np, FORMAT='(20X,I4)'
readf, un, nxg, FORMAT='(20X,I4)'
readf, un, nzg, FORMAT='(20X,I4)'
readf, un, castgrid, castnested, FORMAT='(20X,I1,X,I2)'
readf, un, xfocus, zfocus, FORMAT='(20X,F,X,F)'
readf, un, nxr, FORMAT='(20X,I4)'
readf, un, nzr, FORMAT='(20X,I4)'
readf, un, castradial, castradialpowerx, castradialpowerz , FORMAT='(20X,I1,X,F,X,F)'
readf, un, xupper, FORMAT='(20X,F6)'
readf, un, zupper, FORMAT='(20X,F6)'
readf, un, xphys, FORMAT='(20X,F6)'
readf, un, scrap
readf, un, totnumph, FORMAT='(20X,I8)'
readf, un, theta1, FORMAT='(20X,F4.1)'
readf, un, theta2, FORMAT='(20X,F4.1)'
readf, un, scrap
readf, un, maxscats, FORMAT='(20X,I8)'
readf, un, lymanswitch, FORMAT='(20X,I1)'
readf, un, randomizeangle, FORMAT='(20X,I1)'
readf, un, usedustabs, FORMAT='(20X,I1)'
readf, un, usedustsca, FORMAT='(20X,I1)'
readf, un, constabs, FORMAT='(20X,D)'
readf, un, constsca, FORMAT='(20X,D)'
readf, un, constepsilon, FORMAT='(20X,D)'
readf, un, constg, FORMAT='(20X,D)'
readf, un, usemolabs, FORMAT='(20X,I1)'
readf, un, xrayflag, FORMAT='(20X,I1)'

readf, un, scrap

readf, un, plotglobal, FORMAT='(20X,I1)'
readf, un, plotdiskinit, FORMAT='(20X,I1)'
readf, un, plotrun, FORMAT='(20X,I1)'
readf, un, plotirr, FORMAT='(20X,I1)'
readf, un, plotreg, FORMAT='(20X,I1)'
readf, un, plotorig, FORMAT='(20X,I1)'
;readf, un, usedustsca, FORMAT='(20X,I1)'

close, un

;===============================

file.save=savfile
file.disk=diskfile

para.wavelength=wavelength
para.frequency=const.c/(para.wavelength*1E-8)

para.starlumin=double(starlumin)*1.D44
para.isrfg=isrfg   ; isrfg is in Habings.

grid.npoisson=np
flags.castgrid=castgrid
flags.castnested=castnested
flags.castradial=castradial
flags.castradialpowerx=castradialpowerx ; goes in Flags structure, as value does double-duty as a flag.
flags.castradialpowerz=castradialpowerz

grid.nxg=nxg
grid.nzg=nzg
flags.xf=xfocus  ; goes in Flags structure, as value does double-duty as a flag.
flags.zf=zfocus

grid.nxr=nxr
grid.nzr=nzr

grid.xupper=xupper
;grid.yupper=yupper
grid.zupper=zupper
grid.xphys=xphys
grid.xphyscm=grid.xphys*const.auincm

para.totnumph=totnumph
para.theta1=theta1
para.theta2=theta2

para.maxscats=maxscats
flags.lyman=lymanswitch
flags.randomizeAngle=randomizeangle
flags.useabs=usedustabs
flags.usesca=usedustsca

para.constabs=constabs
para.constsca=constsca
para.consteps=constepsilon
para.constg=constg

flags.usemolabs=usemolabs

;===================  plotting flags.  Master over-ride

plotFlags.global=plotglobal

if(plotFlags.global EQ 1) then begin
  plotFlags.diskInit=plotdiskinit
  plotFlags.overplotmain=plotrun
  if(plotFlags.diskInit EQ 0 AND plotFlags.overplotmain EQ 1) then plotFlags.overplotmain=0
  plotFlags.irr=plotirr
  plotFlags.gridded=plotreg
  plotFlags.origgrid=plotorig
  plotFlags.avsimple=0
endif else begin
  plotFlags.diskInit=0
  plotFlags.overplotmain=0
  plotFlags.irr=0
  plotFlags.gridded=0
  plotFlags.origgrid=0
  plotFlags.avsimple=0
endelse


;===============================

print, '======================='
print, 'Save file: ',file.save
print, 'Input disk file: ',file.disk
print, ''
print, 'wavelength: ',para.wavelength
; print, 'frequency : ',para.frequency
print, 'Grid Data:'
print, 'grid.nxg: ',grid.nxg
print, 'grid.nzg: ',grid.nzg
print, 'flags.castgrid: ',flags.castgrid
print, 'flags.castnested Nlevels: ',flags.castnested
print, 'grid.nxr: ',grid.nxr
print, 'grid.nzr: ',grid.nzr
print, 'flags.castradial, flags.castradialpowerx, flags.castradialpowerz',flags.castradial, flags.castradialpowerx,flags.castradialpowerz
print, ''

print, 'grid.xupper: ',grid.xupper
print, 'grid.zupper: ',grid.zupper
print, ''

print, 'para.totnumph: ',para.totnumph
print, 'para.maxscats: ',para.maxscats
print, 'flags.lyman: ',flags.lyman
print, 'flags.randomizeAngle: ',flags.randomizeAngle
print, 'flags.useabs: ',flags.useabs
print, 'flags.usesca: ',flags.usesca

;print, 'Changing working directory to :',prjdir
;stop
end
