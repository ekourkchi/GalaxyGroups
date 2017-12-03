;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
Function pixel_vertices, ipring, nside

pix2ang_ring, nside, ipring, theta, phi
latlon = [90.-(theta*180./!pi), phi*180./!pi]  ; making astronoomy sense of coordinates

pix2vec_ring, nside, ipring, vector, vertex
latlon = [latlon, vertices(vertex)]

; latlon = [lat_c, lon_c, lat_v1, lon_v2, ...,lat_v4,lon_v4]
; lat_c and lon_c are coordinates of the pixel center
; lat_v1 ... are the coordinates of 4 vertices
return, latlon

END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

Function vertices, vertex

latlon = []

for i=0, 3 do begin
  vec2ang, vertex[0,*,i], lat, lon, /astro
  latlon = [latlon, lat]
  latlon = [latlon, lon]

endfor


return, latlon

END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;




PRO dust_test_v01, galactic=galactic, supergalactic=supergalactic


  
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

gall_v1 = []
galb_v1 = []

gall_v2 = []
galb_v2 = []

gall_v3 = []
galb_v3 = []

gall_v4 = []
galb_v4 = []


nside  = 128L


IF KEYWORD_SET(galactic) THEN BEGIN
  

  for ipring=0, 12*nside*nside-1 do begin

    latlon = pixel_vertices(ipring, nside)

    gall = [gall, latlon[1]]
    galb = [galb, latlon[0]]
    

    
    gall_v1 = [gall_v1, latlon[3]]
    galb_v1 = [galb_v1, latlon[2]]    
    
    gall_v2 = [gall_v2, latlon[5]]
    galb_v2 = [galb_v2, latlon[4]]
    
    gall_v3 = [gall_v3, latlon[7]]
    galb_v3 = [galb_v3, latlon[6]]
    
    gall_v4 = [gall_v4, latlon[9]]
    galb_v4 = [galb_v4, latlon[8]]
    
  endfor  
  
  
  print, n_elements(gall)
;   print, gall


  ipath = '/home/ehsan/PanStarrs/glga/dust/maps/'
  value = wcs_getval( $
	['SFD_dust_4096_ngp.fits', 'SFD_dust_4096_sgp.fits'], $
	gall, galb, path=ipath, interp=interp, noloop=noloop, verbose=verbose)

	
  outstr={GL:0.d,GB:0.d,GL_v1:0.d,GB_v1:0.d,GL_v2:0.d,GB_v2:0.d,GL_v3:0.d,GB_v3:0.d,GL_v4:0.d,GB_v4:0.d,EBV:0.d}

  out=replicate(outstr,n_elements(value))

  out.GL  = gall
  out.GB  = galb
  
  out.GL_v1  = gall_v1
  out.GB_v1  = galb_v1
    
  out.GL_v2  = gall_v2
  out.GB_v2  = galb_v2
  
  out.GL_v3  = gall_v3
  out.GB_v3  = galb_v3
  
  out.GL_v4  = gall_v4
  out.GB_v4  = galb_v4
        
  out.EBV = value

  mwrfits,out,'EBV.nside.128.gal.fits',/create
ENDIF

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;








  
END









