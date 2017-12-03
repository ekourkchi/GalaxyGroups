pro aitoff_example,hardcopy=hardcopy

linethick=5
charthick=7
device,retain=2

if (keyword_set(hardcopy)) then PS_start,filename='fig.ps',xsize=32/2.54 , ysize=19.0/2.54, /inches
theta=findgen(30)/29*2*!pi 
usersym,sin(theta),cos(theta),/fill


;-------------------------------------------plotting the declination signs and frame -------------------------------------
centerl=0.;in degrees
centerb=0.;in degrees
chipwidth=2/3.;in degrees
chiplength=2/3.;in degrees

 leftsidex=[centerl-chipwidth/2,centerl-chipwidth/2]
 leftsidey=[centerb-chiplength/2,centerb+chiplength/2]
 
 rightsidex=[centerl+chipwidth/2,centerl+chipwidth/2]
 rightsidey=[centerb-chiplength/2,centerb+chiplength/2]
 
 topsidex=[centerl-chipwidth/2,centerl+chipwidth/2]
 topsidey=[centerb+chiplength/2,centerb+chiplength/2]
 
 bottomsidex=[centerl-chipwidth/2,centerl+chipwidth/2]
 bottomsidey=[centerb-chiplength/2,centerb-chiplength/2]

aitoff,leftsidex,leftsidey,x,y
;-----------------------------------------------------------------------------------------------

readcol,'points.dat',ra_s1,dec_s1
xs1=ra_s1+90.
ys1=dec_s1

aitoff,xs1,ys1,xall,yall
PLOTSYM,0,thick=4.
oplot,xall,yall,psym=8,color=cgcolor(' orange'),symsize = 1.9



;------------------------------------------- the signs -------------------------------------

 aitoff,180.,0,xt,yt
 xt=-xt
 xyouts,xt-15,yt-2,'18h',charsize=1.4,charthick=1.9  ; 300 deg
 xyouts,xt+1,yt-2,' ',charsize=1.4,charthick=1.9     ; 
 xyouts,xt+180.,yt-5.,'6h',charsize=1.4,charthick=1.9 ;120 deg
aitoff,0,30.,xt,yt
xt=-xt
xyouts,xt,yt,'$30^{\circ}$',charsize=1.4,charthick=1.9
aitoff,0,-30.,xt,yt
xt=-xt
xyouts,xt,yt,'$-30^{\circ}$',charsize=1.4,charthick=1.9
aitoff,0,60.,xt,yt
xt=-xt
xyouts,xt,yt,' 60^{\circ}',charsize=1.4,charthick=1.9
aitoff,0,-60.,xt,yt
xt=-xt
xyouts,xt,yt,'-60^{\circ}',charsize=1.4,charthick=1.9
aitoff_grid



if (keyword_set(hardcopy)) then PS_end,/png

stop
end
