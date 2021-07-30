pro uvrt2dcyl, INPUT_FILE= paraFile, NO_DEVICE=nodevice

;==========================================================
; Radiative transfer in 2D Cylindrical. version 1.3
; thomas Bethell. Apr 2010
;
;v1.1
;  put all global data_structures into a common block.  
;  Common blocks are how global variables are communicated in IDL.
;v1.2
;  identity of parameter file is specified at runtime.
;  added radial grid points.
;  made some of the loops LONG variables.
;v1.3
;  added NO_DEVICE keyword.  More powerful device suppression than plotFlags.Global flag.
;  Setting this allows code to run remotely without any call to things like windows, or device properties.
;v1.4
;  Seem to have fixed dying photon problem.  Added ISRF option -> New line in parameter file.


common GlobalVars, plotFlags, const, flags, para,$
      file, cell, grid, vert, tr, ph, diag, out, orig,$
      levs, xrayflag

;===========================================================
;  Windows
;===========================================================

 ; window,0
;	set_plot, 'x'
	;!p.multi=[0,1,1]
; xrayflag = 1
;===========================================================
;  Set working directory, fix display to retain display 
;  image in proper colors
;===========================================================

 ; cd, '/n/Users/tbethell/IDLWorkspace/rt2dcyl/'  ; invoke this once in IDLDE if needed

;print,'dsds' 
;stop 
 
 
;===========================================================
;  Flags
;===========================================================

	flags={castgrid:0,$
	castnested:0,$
	castradial:0,$
	castradialpowerx:0.0,$
	castradialpowerz:0.0,$
	xf:0.0,$
	zf:0.0,$
	lyman:0,$
	useabs:0,$
	usesca:0,$
	usemolabs:0,$
	randomizeAngle:0}
	
	plotFlags={nodevice:0,$
	global:1,$
	diskinit:1,$
	overplotmain:1,$    ; 2=no overplot, 1=overplotmain, 3=plot individual cells?
	irr:0,$ 
  gridded:0,$
  origgrid:0,$
	avsimple:0}
  
  if(KEYWORD_SET(nodevice)) then plotFlags.nodevice=1
  
  levs={numl:101,$
        numlfew:11,$
        nh:10.^(8.5*findgen(101)/100.),$
        nhfew:10.^(8.5*findgen(11)/10.),$
        j:10.^(8.5*findgen(101)/100.),$
        nph:10.^(8.5*findgen(101)/100.),$
        nphfew:10.^(8.5*findgen(11)/10.),$
        a:findgen(101)/100.,$
        f:10.^(8+8*findgen(101)/100.)}

  diag={fup:0L}     ; diagnostics

;===========================================================
;  files
;===========================================================

  file={save:"k",disk:"d"}

;===========================================================
;  run parameters
;===========================================================

para={wavelength:0.0,frequency:0.0D,starlumin:0.0D,isrfg:0.0D,maxscats:0L,totnumph:0L,$
        theta1:0.0,theta2:0.0,minwabs:0.0,$
        constabs:0.0,constsca:0.0,constg:0.0,consteps:0.0}

;===========================================================
;  time-keeping etc.
;===========================================================

time={veryStart:0,$
      start:0,$
      veryEnd:0}

time.veryStart=systime(1, /seconds)

;===========================================================
;  basic physical constants etc.
;===========================================================

const={c:2.997924E10,$   ;[cm s^-1]
      rootpi:1.77245385,$
      esq:2.307E-19,$
      me:9.109E-28,$
      mh:1.67E-24,$
      kb:1.3807E-16,$
      GravConst:6.673E-8,$
      Msol:3.9E33,$
      auincm:double(1.495980000000D13),$
      degtorad:2.0D*!dpi/360.0D	,$
      freql:9.936E7,$
      angcen:1215.67 ,$   ; wavelength of lyman alpha Angs
      freqcen:2.466E15 ,$     ;2.46607E15 
      lyaenergy:10.4D,$  ; Lya E in eV
      dopunit:0,$
      f:0.4162,$
      ftrue:0.4162,$
      voigtatomic:0.0,$
      ergtoev:624150974451.2D,$
      lsolarerg:3.9D33,$
      lsolarev:0}   ; ergs s-1

	;const.dopunit=(1.E5/const.c)*const.freqcen   ; unit doppler is for 1km s-1.  Units Hz
  const.freqcen=const.c/(const.angcen*1.E-8)      ;2.466E15
  
  const.voigtatomic=const.f*const.rootpi*const.esq/(const.me*const.c)
  
	const.lsolarev=const.lsolarerg*const.ergtoev
	;const.lsolarlyan=const.lsolarev/const.lyaenergy   ; Solar equivalent lyalpha photons s^-1
  
  ;===========================================================
  ;  grid
  ;===========================================================
  
  grid={nxg:0.0,nzg:0.0,nxr:0.0,nzr:0.0,xphys:0.0,xphyscm:0.0,xgap:0.0,xupper:0.0,zupper:0.0,$
        px:0.0,pz:0.0,npoisson:0L}
  
  ;===========================================================
  ;  get device properties
  ;===========================================================
 
 if (plotFlags.nodevice EQ 0) then begin
  Device, Get_Visual_Name=theVisual, Get_Visual_Depth=theDepth
  Print, 'Setting window display device ... ',theVisual, theDepth
  
  device, retain=2
  device, decomposed=0
 endif else begin
 
 Print, 'No devices will be used.  No plots to eps or windows.'
 
 endelse
 
 ;stop
 
  ;===========================================================
  ;  Read Basic parameters
  ;===========================================================

if KEYWORD_SET(paraFile) then begin
  a_readparas, paraFile
endif else begin
  a_readparas, "parameters.txt"
endelse

;a_readparas, paraFile

	;===========================================================
	;  Copy parameter file.  Add snippet of disk file.  Add textual summary
	;===========================================================

;if KEYWORD_SET(paraFile) then begin
;  a_copyparas, paraFile
;endif else begin
;  a_copyparas, "parameters.txt"
;endelse

	;===========================================================
	; Load in extinction curve.  
	; Either construct one, load premade or use Nuria's two component dust proprties
	;===========================================================

	;window, 1, xs=800,ys=400
	;wset, 1
	;!p.multi=[0,2,1]

	;print, 'Loading dust extinction curve ...'

	;===========================================================
	; Load in water photodissociation cross-section
	;===========================================================

	;window, 1, xs=800,ys=400
	;wset, 1
	;!p.multi=[0,2,1]

	;print, 'Loading water UV cross-sections...'
	;A_load_water_crossx

		;=============================================================
		;  Make grid.  Fill it with density, dust to gas, h2frac,$
		;  cell vols/areas,
		;=============================================================

		print, 'Loading disk and optional background grid ...'

	
		a_disk_create, xp, zp, nhp, ntr, comx, comz, nh, epsi, sca, abss, g, fh, tg, fwater, foh, trang, trangx, trangz, trar, trvol
  
    ;===========================================================
    ; just as for generator nodes, check for comx,comz inside gap
    ; note: code units.
    ;===========================================================
      
    whingap=where(comx LE grid.xgap, numgap)
    if(numgap gt 0) then  nh[whingap]=0.0
 
 
 ;======== look at innermost cells for nonzero error
 
 ;print, nh[whingap]
 ;stop
 
    ;===========================================================
    ; Now disk vertices are known, define vert, cell, and trang block
    ;===========================================================
    
    vert={x:xp,z:zp,nh:nhp}
    
    cell={ntr:ntr,$
      comx:comx,$
      comz:comz,$
      nh:nh,$
      tg:tg,$
      td:nh*0.,$
      fh:fh,$
      fwater:fwater,$
      foh:foh,$
      abs:abss,$
      sca:sca,$
      epsi:epsi,$
      g:nh*0.+0.5,$
      dop:nh*0.}
 
 ;======== DEBUG
 
; whhlayer=where(cell.nh LT 1E5 AND cell.nh GT 1E4)
;cell.fh[whhlayer]=1.0
;cell.tg[*]=1000.0
 
      
;      cell.abs[*]=1E-21
;      cell.sca[*]=1E-21
;      cell.epsi[*]=1
;      
     ph={nph:DOUBLE(nh*0.0),$
         j:nh*0.,$
         l:nh*0.,$
         n:nh*0.,$
         fx:nh*0.,$
         fz:nh*0.,$
         nphnorm:nh*0.,$
         jnorm:nh*0.,$
         fxnorm:nh*0.,$
         fznorm:nh*0.,$
         a:nh*0.,$
         packetlumin:0.0D}

    tr={ntr:ntr,$
          arr:trang,$
          x:trangx,$
          z:trangz,$
          ar:trar,$
          vol:trvol,$
          neigh:trang*0.,$
          edgetype:trang*0.}
 
	  ;=================================================
    ;  plot
    ;=================================================

 if(plotFlags.nodevice EQ 0 AND plotFlags.global EQ 1 AND plotFlags.diskInit EQ 1) then begin
      print, 'Making a plot of the disk'
      a_plotquad
    endif

;stop

		;=================================================
		;  Finding neighbours
		;=================================================

		print, 'Finding neighbours by using Delaunay Tesselation...'

  a_neighbours

	;===========================================================
	; Given numtotph and angles, compute lumin per packet.
	;===========================================================

	a_compute_packet_luminosity

	;===========================================================


;
;
;
;
;
;
;
;
;
;
;
;
;
;
;
;

;plotFlags.overplotmain=1

time.start=systime(1, /seconds)

	;=================================================
	;  ffor various values, if selected in parameter file
	;=================================================

	;	A_wav_interp, 0.1215, exwav, exabs, exsca, exg, absdust, scadust, g

   if(para.constabs NE 0) then cell.abs[*]=para.constabs
   if(para.constsca NE 0) then cell.sca[*]=para.constsca
   if(para.consteps NE 0) then cell.epsi[*]=para.consteps
   if(para.constg NE 0) then cell.g[*]=para.constg

	 if(flags.useabs eq 0)  then cell.abs[*]=0.
	 if(flags.usesca eq 0)  then cell.sca[*]=0.
	  ;if(flags.usemolabs eq 0)  then cell.molabs[*]=0.
  	
;print, 'debug'
;print, cell.abs[1],cell.sca[1],cell.g[1],cell.epsi[1]
;stop

		;=================================================
		; and watercrossx at 1216A
		;=================================================
	
		;water.crossxatwav=water.crossx(max(where(water.wavm le 0.1216)))

		;=================================================
		;  
		;=================================================
				
		;const.f=const.ftrue

		;if(flags.lyman ne 1) then const.f=0.   ; merely turns off resonant scattering, 
						;not Ly-a luminosity (that is what dolyman flag is for).

	;=================================================
	; Set up arrays.  Doppler unit in Hz
	;=================================================

	cell.dop=0.*cell.nh+((sqrt(2.*const.kb*cell.tg/const.mh))/const.c)*const.freqcen   ; Hz. Delta nu sub D.
	 
 	;=================================================
	;  Run main2 which does RT.
	;=================================================

	main2

	;==========================================================
	;  Photon density.  n per cm^3 scaled for Lstar(4pi sterad)=1
	;==========================================================

	a_photon_density
	
;	endif

time.veryEnd=systime(1,/seconds)	

;===========================================
; Vertical (or otherwise) A_V for every cell 
;===========================================

;	A_visual_extinction, scadustv, absdustv
	
;	A_avcell, scadustv, absdustv
 ; cell.av=0.*dcell

;==========================================================
; print out some checks
;==========================================================

	;a_print_some_values

;=============================================
;  Info  run time etc.
;=============================================

  totaltime=(time.veryEnd-time.start)

  print, 'Raytracing runtime (secs, hours):',totaltime,totaltime/3600.
  ;print, 'Total runtime (secs)     :',endtime-verystarttime

;==========================================================
;  Recast full results onto original grid.
;==========================================================
    
  a_recast_onto_orig
  a_save_results

;==========================================================
;  Plot, regrid results.
;==========================================================

  MAIN3

;=============================================
;  save important results in txt file, matching the input disk file.
;=============================================

;print, 'Total runtime (secs, over all iterations) :',systime(1,/seconds)-time.veryStart
print, 'Finished.'

end
