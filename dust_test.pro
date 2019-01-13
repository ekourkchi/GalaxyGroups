PRO dust_test, galactic=galactic, supergalactic=supergalactic


; xrange1 = [0, 360.]
; yrange1 = [-30., 30.]
; 
; 
; p = PLOT(/nodata,xrange1,yrange1, thick=5,xrange=xrange1,yrange=yrange1,$
; 	xtitle='Gl',ytitle='Gb',$
; 	xstyle=1,ystyle=1,title='Galactic Box',xthick=2,ythick=2, DIMENSIONS=[1800, 300], FONT_SIZE=16)
; 
; ctable = COLORTABLE(10)
; 
; delta = 3.

; for gl = xrange1[0], xrange1[1], delta do begin
; for gb = yrange1[0], yrange1[1], delta do begin
; 
; 
; 
; value = dust_getval(gl+0.5*delta, gb+0.5*delta)
; ind = floor(255-value*100)
; if ind lt 0 then ind=0
; if ind gt 255 then ind=255
; 
; poly = POLYGON([gl,gl,gl+delta,gl+delta,gl],[gb,gb+delta,gb+delta,gb,gb], $
;    /DATA, /FILL_BACKGROUND, $
;   FILL_COLOR=ctable[ind], LINESTYLE=6)
; 
;   
; print, gl, gb, value, ind
; 
;   endfor
;   endfor

  
gall = []
galb = []
sgall = []
sgalb = []
deg = 0.5
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


IF KEYWORD_SET(supergalactic) THEN BEGIN

  for sgl=0., 360.-deg, deg do begin
  for sgb=-90., 90.-deg, deg do begin
    glactc, ra,dec,2000,sgl,sgb,2,/SuperGalactic
    glactc, ra,dec,2000,gl,gb,1
    
    gall = [gall, gl]
    galb = [galb, gb]
    
    sgall = [sgall, sgl]
    sgalb = [sgalb, sgb]
    
  endfor
  endfor


  ipath = '/home/ehsan/PanStarrs/glga/dust/maps/'
  value = wcs_getval( $
	['SFD_dust_4096_ngp.fits', 'SFD_dust_4096_sgp.fits'], $
	gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)

	
  outstr={SGL:0.d,SGB:0.d,EBV:0.d}

  out=replicate(outstr,n_elements(value))

  out.SGL  = sgall
  out.SGB  = sgalb
  out.EBV = value

  mwrfits,out,'EBV.0.8.deg.fits',/create
ENDIF 


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

gall = []
galb = []
sgall = []
sgalb = []
deg = 0.5


IF KEYWORD_SET(galactic) THEN BEGIN
  for gl=0., 360.-deg, deg do begin
  for gb=-90., 90.-deg, deg do begin
    
    gall = [gall, gl]
    galb = [galb, gb]
    
  endfor
  endfor


  ipath = '/home/ehsan/PanStarrs/glga/dust/maps/'
  value = wcs_getval( $
	['SFD_dust_4096_ngp.fits', 'SFD_dust_4096_sgp.fits'], $
	gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)

	
  outstr={GL:0.d,GB:0.d,EBV:0.d}

  out=replicate(outstr,n_elements(value))

  out.GL  = gall
  out.GB  = galb
  out.EBV = value

  mwrfits,out,'EBV.0.5.gal.fits',/create
ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  
END

