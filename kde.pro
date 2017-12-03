;+
; NAME:
;    KDE
;
; PURPOSE:
;    Estimate the probability density underlying a set of discrete
;    samples (measurements) using the kernel density estimator method.
;
; CATEGORY:
;    Statistics
;
; CALLING SEQUENCE:
;    d = kde(x, t)
;
; INPUTS:
;    x: discrete samples of the desired distribution.
;
;    t: values for which the probability density is required.
;
; KEYWORD PARAMETERS:
;    scale: smoothing factor, also called the bandwidth
;        used to compute the kernel density estimate
;
;    weight: weighting for sampled points.
;
; KEYWORD FLAGS:
;    By default, KDE uses the Epanechnikov kernel to compute
;    the kernel density estimate, this can be overridden by
;    setting one of the following flags:
;    GAUSSIAN: use Gaussian kernel
;    TRIANGULAR: use triangular kernel
;    BIWEIGHT: use biweight kernel
;
; OUTPUTS:
;    d: probability density estimated at each value specified by t
;
; RESTRICTIONS:
;    Gaussian kernel used for multivariate systems.
;
; PROCEDURE:
;    d_i = (1/n) \sum_{j = 1}^n K((x_j - t_i)/h) / h
;
;    where h is the estimated optimal smoothing parameter and
;    where K(z) is the selected kernel.
;
; REFERENCE:
; B. W. Silverman,
; Density Estimation for Statistics and Data Analysis
; (CRC Press, Boca Raton, 1998)
;
; EXAMPLE:
;    IDL> x = randomn(seed, 1000)
;    IDL> t = 2. * findgen(100)/99.
;    IDL> d = kde(x,t)
;    IDL> plot, t, d
;    IDL> plot, t, histogram(x, min=0, max=2, nbins=100), /noerase
;
; MODIFICATION HISTORY:
; 09/18/2010 Written by David G. Grier, New York University
; 10/08/2011 DGG: Corrected for case when more than 3/4 of input data
;     have same value.  Thanks to Dan Hartung (UW Madison) 
;     for bringing this bug to light.
; 10/18/2011 DGG: Added WEIGHT keyword.
;
; Copyright (c) 2010, 2011 David G. Grier
;
;
; UPDATES:
;    The most recent version of this program may be obtained from
;    http://physics.nyu.edu/grierlab/software.html
; 
; LICENSE:
;    This program is free software; you can redistribute it and/or
;    modify it under the terms of the GNU General Public License as
;    published by the Free Software Foundation; either version 2 of the
;    License, or (at your option) any later version.
;
;    This program is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
;    General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with this program; if not, write to the Free Software
;    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
;    02111-1307 USA
;
;    If the Internet and WWW are still functional when you are using
;    this, you should be able to access the GPL here: 
;    http://www.gnu.org/copyleft/gpl.html
;-

function kde_nd, x, y, $
                 weight = weight, $
                 scale = scale

COMPILE_OPT IDL2, HIDDEN

sx = size(x, /dimensions)
sy = size(y, /dimensions)

nd = sx[0]                      ; number of dimensions
nx = sx[1]                      ; number of data points
ny = sy[1]                      ; number of sampling points

if n_elements(weight) ne nx then $
   weight = replicate(1., nx)

; optimal smoothing parameter in each dimension
; Silverman Eqs. (3.30) and (3.31)
sx = stddev(x, dimension=2)
rx = fltarr(nd)
for d = 0, nd-1 do $
   rx[d] = iqr(x[d,*])
h = sx
w = where(rx gt 1e-10, ngood)
if ngood gt 0 then $
   h[w] = h[w] < rx[w]/1.34
h *=  0.9 / nx^0.2

scale = h

; density estimate
; Silverman Eq. (2.15) and Table 3.1
fac = replicate(1., nx)
res = fltarr(ny)
hfac = h # fac

for j = 0, ny-1 do begin
   z = total(0.5 * ((x - y[*,j] # fac) / hfac)^2, 1)
   w = where(z lt 20, ngood)
   if ngood gt 0 then $
      res[j] = total(weight[w]*exp(-z[w])) / nx
endfor

res /= 2. * !pi * total(h^2)^(nd/2) ; normalization

return, res
end

function kde_1d, x, y, $
                 weight = weight, $
                 scale = scale, $
                 biweight = biweight, $
                 triangular = triangular, $
                 gaussian = gaussian

COMPILE_OPT IDL2, HIDDEN

nx = n_elements(x)              ; number of data points
ny = n_elements(y)              ; number of samples

if n_elements(weight) ne nx then $
   weight = replicate(1., nx)

; optimal smoothing parameter
; Silverman Eqs. (3.30) and (3.31)
sx = stddev(x)                  ; standard deviation
rx = iqr(x)                     ; interquartile range
if rx lt 1e-10 then $           ; more than 3/4 data have same value
   h = 0.9 * sx / nx^0.2 $
else $
   h = 0.9 * (sx < rx/1.34) / nx^0.2
scale = h

; density estimate
; Silverman Eq. (2.15) and Table 3.1
t = x/h
s = y/h
res = fltarr(ny)                ; result

if keyword_set(biweight) then begin
   for j = 0, ny-1 do begin
      z = (t - s[j])^2
      w = where(z lt 1., ngood)
      if ngood gt 0. then $
         res[j] = 15.*total(weight[w]*(1.-z[w])^2)/(16.*h*nx)
   endfor
endif $                     
else if keyword_set(triangular) then begin
   for j = 0, ny-1 do begin
      z = abs(t - s[j])
      w = where(z lt 1., ngood)
      if ngood gt 0 then $
         res[j] = total(weight[w]*(1. - z[w]))/(h*nx)
   endfor
endif $                     
else if keyword_set(gaussian) then begin
   for j = 0, ny-1 do begin
      z = 0.5 * (t - s[j])^2
      w = where(z lt 20, ngood)
      if ngood gt 0 then $
         res[j] = total(weight[w]*exp(-z[w]))/(sqrt(2.*!pi)*h*nx) 
   endfor
endif $
else begin                      ; Epanechnikov
   for j = 0, ny-1 do begin
      z = (t - s[j])^2
      w = where(z lt 5, ngood)
      if ngood gt 0 then $
         res[j] = 0.75*total(weight[w]*(1.-z[w]/5.))/(sqrt(5.)*h*nx) 
   endfor
endelse

return, res
end

function kde, x, y, $
              weight = weight, $
              scale = scale, $
              gaussian = gaussian, $
              biweight = biweight, $
              triangular = triangular

COMPILE_OPT IDL2
  
sx = size(x)
sy = size(y)

if sx[0] gt 2 then $
   message, "data must be organized as [ndims,npoints]"

if sy[0] ne sx[0] then $
   message, "inputs must have the same number of dimensions"

if (sx[0] eq 2) and (sx[1] ne sy[1]) then $
   message, "inputs must have the same number of dimensions"

ndims = (sx[0] eq 2) ? sx[1] : 1

if ndims gt 1 then begin
   if keyword_set(biweight) or keyword_set(triangular) then $
      message, "Multidimensional: using Gaussian kernel", /inf
   res = kde_nd(x, y, weight = weight, scale = scale)
endif else $
   res = kde_1d(x, y, weight = weight, scale = scale, $
               gaussian = gaussian, $
               biweight = biweight, $
               triangular = triangular)

return, res
end


