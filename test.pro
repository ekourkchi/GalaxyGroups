pro test

a=mrdfits('fsat_60000.fits',1)
fsat=a.fsat
;readcol, 'dpeak_10000.dat', fsat
fsat=(fsat+1.5e-5)*15

step=0.1e-4
min=0.5e-4
range=0.001
bins=CEIL(range/step)
fsatbin=fltarr(bins)
for i=0,(bins-1) do begin
    fsatbin[i]=min+(step/2)
    min=min+step
endfor

PDF = kde(fsat,fsatbin, /gaussian)

PDF = PDF/total(PDF*step)
plot, fsatbin*10000.,pdf,xrange=[2,12]

 
for i=0,10 do begin

print, fsatbin[i]*10000., pdf[i]

endfor


g = gaussfit(fsatbin[0:100]*10000., pdf[0:100]  ,a, NTERMS=3 )
print,a


x=0
y=0

z=2
while z le 12 do begin


x=[x,z]
u = (z-a[1])/a[2]
y=[y, a[0]*exp(-u*u/2.)]


z+=0.1
endwhile


;oplot, x,y



end