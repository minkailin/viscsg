FUNCTION logticks_exp, axis, index, value
   ; Determine the base-10 exponent
   exponent   = LONG( ALOG10( value ) )
   ; Construct the tickmark string based on the exponent
;   tickmark = '10!E' + STRTRIM( STRING( exponent ), 2 ) + '!N'
    tickmark = textoidl('10^{' + STRTRIM( STRING( exponent ), 2 )+'}')
   ; Return the formatted tickmark string
   RETURN, tickmark
END

pro result, xrange=xrange, yrange=yrange, label=label, legend=legend, title=title , ytickinterval=ytickinterval, yminor=yminor

!p.font = 0 

params = dblarr(7,1)
openr,1,'params.dat'
readf,1,params
close,1
bigQ   = params(0,0)
gmma   = params(1,0)
mu     = params(2,0)
lambda = params(3,0)
smallq = params(4,0)
bulkv  = params(5,0) 
tirr   = params(6,0) 

lines = file_lines('output.dat')
array = dblarr(7,lines)

openr,1,'output.dat'
readf,1,array
close,1

baxis = array(0,*)
alpha = array(1,*)
rate  = array(2,*)
freq  = array(3,*)
kmax  = array(4,*)
error = array(5,*)
Qarr  = array(6,*)


temp = min(abs(alpha-1d0),agrid)

temp = min(abs(1d0/(2d0*!dpi*10.)-rate),rgrid) 

;special cases with alpha-beta relation 
;2D only 

rate_special = rate
kmax_special = rate 

mink = 1d-1
maxk = 1d1
nk   = 8192
kaxis = mink + (maxk-mink)*dindgen(nk)/(nk-1d0) 
rate_special_fixb = dindgen(nk)
for i=0, lines-1 do begin

;for mu=-1, lambda=1 (standard alpha disk) 
;a4 = 1d0
;a3 = (7d0/3d0)*alpha(i)*k^2
;a2 = 2d0*(2d0-smallq)+gmma*k^2 - 2d0*k/bigQ + (4d0/3d0)*alpha(i)^2*k^4 + 2d0*(alpha(i)*smallq*k)^2*(gmma-1d0)
;a1 = gmma*k^2 - 2d0*k/bigQ + 2d0*gmma*smallq + (gmma-1d0)*( 8d0*(alpha(i)*smallq*k)^2/3d0 - 2d0*smallq*(2d0-smallq) )
;a1*= alpha(i)*k^2
;a0 = -4d0*(alpha(i)*smallq*k)^2*(gmma-1d0)*k/bigQ 

;for mu=-1, lambda=0 (gammie's disk)
;; a0 = -((2.*alpha(i)^2*(-1. + gmma)*k^2*smallq^2*k)/bigQ)
;; a1 = (alpha(i)*(-1. + $ 
;;                 gmma)*smallq^2*(bigQ*(k^4*(4.*alpha(i)^2 + (3.*gmma)/((-1. + gmma)*smallq^2)) - $
;;                                       6.*(-2. + smallq) + (6.*k^2*(-2. + smallq))/smallq) - $
;;                                 6.*(1. + k^2/((-1. + gmma)*smallq^2))*k))/(3.*bigQ)
;; a2 = (bigQ*(12. + 4.*alpha(i)^2*k^4 - 6.*smallq - 7.*alpha(i)^2*k^2*smallq^2 + $
;;             gmma*k^2*(3. + 7.*alpha(i)^2*smallq^2)) - 6.*k)/(3.*bigQ)
;; a3 = alpha(i)*(-1. + gmma)*(1. + (7.*k^2)/(3.*(-1. + gmma)*smallq^2))*smallq^2
;; a4 = 1d0

;max growth rate, inviscid
;a0 = - 1d0/gmma/bigQ^2/baxis(i) 
;a1 = 2d0*(2d0-smallq) - 1d0/gmma/bigQ^2
;a2 = 0d0 
;a3 = 1d0 
;roots = fz_roots([a0,a1,a2,a3]) 
;rate_special(i) = max(real_part(roots))

;estimate for max rate for gammie's disk (mu = -1, lambda=0)
bvisc = bulkv*alpha(i)
kopt = 2d0*(2d0 - smallq) 
kopt/= alpha(i)*(4d0*alpha(i)/3d0 + bvisc + gmma*baxis(i))
kopt = kopt^(1d0/4d0)  

kmax_special(i) = kopt 


rate_special(i) = 2d0/Qarr(i)/kopt - tirr 
rate_special(i)/= 4d0*alpha(i)/3d0 + bvisc + gmma*baxis(i) 

for j=0, nk-1 do begin
k = kaxis(j)
top = 2d0*k/Qarr(i) - gmma*k*k + 2d0*(2d0-smallq)*smallq*(gmma-1d) 
top*= alpha(i)*k*k
bot = gmma*k*k + 2d0*(2d0-smallq) - 2d0*k/Qarr(i) 
rate_special_fixb(j) = top/bot
endfor

;rate_special(i) = max(rate_special_fixb)


endfor 

xtitle = tex2idl('$\beta$') + '!X'
;xtitle = tex2idl('$\alpha$') + '!X'
ytitle = tex2idl('max$(s/\Omega)$') + '!X'

if not keyword_set(title) then title = ''

loadct, 6,/silent
color_arr = dindgen(2)*256d0/2.

set_plot,'ps'
file = strcompress('result.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, baxis,rate,xmargin=[8.3,1.7],ymargin=[3.2,2], ystyle=1, xstyle=1 $
      ,charsize=2, thick=4, xrange=xrange, title=textoidl(title), xtitle=xtitle, yrange=yrange,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog,/ylog, color=color_arr(0), xtickformat='logticks_exp'
oplot, baxis, rate_special, thick=4, color=color_arr(1) 
;oplot, [min(baxis),max(baxis)], [0,0], linestyle=2, thick=2
oplot, [1,1]*baxis(agrid), [1d-8,1d8], thick=2, linestyle=2
xyouts, baxis(agrid)*0.9, 3d-5, textoidl('\alpha>1'), charsize=2, align=1

oplot, [1,1]*baxis[rgrid], [1d-8,1d8], thick=2, linestyle=3
xyouts,baxis[rgrid]*1.2,3d-5, textoidl('t_{grow}>10 orbits'), charsize=2

if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)

    for j=0, n_elements(label)-1 do begin
        ynew = 10.^(alog10(y0) - dy*j)
        oplot, [x0,x1], [ynew,ynew], thick=4, color=color_arr(j) 
      ;  print, ynew
      xyouts, x1, ynew, textoidl(label(j)),charsize=2
   endfor

endif
device,/close

ytitle = tex2idl('$kH$ for max. growth') + '!X'

set_plot,'ps'
file = strcompress('result_kmax.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, baxis,kmax,xmargin=[8.3,1.7],ymargin=[3.2,2], ystyle=0, xstyle=1 $
      ,charsize=2, thick=4, xrange=xrange, title=textoidl(title), xtitle=xtitle, yrange=yrange, yminor=yminor, xtickformat='logticks_exp',  $
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog,color=color_arr(0);,/ylog
oplot, baxis, kmax_special, thick=4, color=color_arr(1)

oplot, [1,1]*baxis(agrid), [1d-8,1d8], thick=2, linestyle=2
xyouts, baxis(agrid)*0.9, 0.2, textoidl('\alpha>1'), charsize=2, align=1

oplot, [1,1]*baxis[rgrid], [1d-8,1d8], thick=2, linestyle=3
xyouts,baxis[rgrid]*1.2,0.2, textoidl('t_{grow}>10 orbits'), charsize=2

if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
   for j=0, n_elements(label)-1 do begin

  ynew = y0 - dy*j
        oplot, [x0,x1], [ynew,ynew], thick=4, color=color_arr(j)
      ;  print, ynew
      xyouts, x1, ynew, textoidl(label(j)),charsize=2
     
   
     endfor
endif

device,/close

ytitle = tex2idl('$\alpha$') + '!X'

set_plot,'ps'
file = strcompress('result_alpha.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, baxis,alpha,xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
      ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=xtitle, yrange=yrange,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog
if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
   for j=0, n_elements(label)-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
endif
device,/close

ytitle = textoidl('1-\alpha\betaq^2(\gamma-1)') + '!X'

Hext = 1d0 - alpha*baxis*smallq^2*(gmma-1d0)

set_plot,'ps'
file = strcompress('result_Hext.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, baxis,Hext,xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
      ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=xtitle, yrange=yrange,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog
if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
   for j=0, n_elements(label)-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
endif
device,/close

ytitle = textoidl('Error') + '!X'

set_plot,'ps'
file = strcompress('result_err.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, baxis,error,xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
      ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=xtitle, yrange=yrange,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog
if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
   for j=0, n_elements(label)-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
endif
device,/close



end
