FUNCTION logticks_exp, axis, index, value
   ; Determine the base-10 exponent
   exponent   = LONG( ALOG10( value ) )
   ; Construct the tickmark string based on the exponent
;   tickmark = '10!E' + STRTRIM( STRING( exponent ), 2 ) + '!N'
    tickmark = textoidl('10^{' + STRTRIM( STRING( exponent ), 2 )+'}')
   ; Return the formatted tickmark string
   RETURN, tickmark
END

pro result, xrange=xrange, yrange=yrange, label=label, legend=legend, title=title 

!p.font = 0 

params = dblarr(7,1)
openr,1,'params.dat'
readf,1,params
close,1
fixQ   = params(0,0)
gmma   = params(1,0)
mu     = params(2,0)
lambda = params(3,0)
smallq = params(4,0)
bulkv  = params(5,0) 
tirr   = params(6,0) 

lines = file_lines('output.dat')
array = dblarr(11,lines)

openr,1,'output.dat'
readf,1,array
close,1

raxis = array(0,*)
alpha = array(1,*)
rate  = array(2,*)
freq  = array(3,*)
kmax  = array(4,*)
error = array(5,*)
Qarr  = array(6,*)
theta = array(7,*)
bcool = array(8,*)
surfd = array(9,*)
temp  = array(10,*) 

;write(10,fmt='(11(e22.15,x))'), raxis(i), alpha_arr(i), max_rate(i), max_freq(i), max_k(i), error(i), Q_arr(i), theta_arr(i), bcool_arr(i), surf_arr(i), temp_arr(i)


temp = min(abs(alpha-0.1d0),agrid)
temp = min(abs(bcool-10d0),bgrid)
print, 'bcool at alpha=0.1 is', bcool(agrid)


temp = min(abs(1d0/(2d0*!dpi*10.)-rate),rgrid) 

;special cases with alpha-beta relation 
;2D only 

rate_special = rate
kmax_special = rate 
for i=0, lines-1 do begin
k = kmax(i)
;estimate for max rate for gammie's disk (mu = -1, lambda=0)

;bvisc = bulkv*alpha(i)
;kopt = 2d0*(2d0 - smallq) 
;kopt/= alpha(i)*(4d0*alpha(i)/3d0 + bvisc + gmma*bcool(i))
;kopt = kopt^(1d0/4d0)  

;kmax_special(i) = 3d0/(2d0*theta(i)*Qarr(i))

;rate_special(i) = 2d0/Qarr(i)/kopt - theta(i)  
;rate_special(i)/= 4d0*alpha(i)/3d0 + bvisc + gmma*bcool(i) 

kmax_special(i) = 3d0/(2d0*theta(i)*Qarr(i))
rate_special(i) = 27d0*alpha(i)/(16d0*theta(i)^3*Qarr(i)^4) 

endfor 




set_plot,'ps'
file = strcompress('result_basic.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=6,xoffset=0,yoffset=0,/inches,/color
multiplot,/reset
multiplot,[1,4],margins=[0.12,0.1,0.05,0.05],rowspacing=0.1
plot, raxis,Qarr,xmargin=[8.3,1.7],ymargin=[5,0], ystyle=1, xstyle=1 $
      ,charsize=1.5, thick=4,yrange=[1,10],$
      linestyle = 0, ytitle = 'Q', xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog;,/ylog, ytickformat='logticks_exp'
multiplot
plot, raxis,alpha,xmargin=[8.3,1.7], ystyle=2, xstyle=1 $
      ,charsize=1.5, thick=4,yrange=[1d-3,1]*0$
      ,linestyle = 0, ytitle = textoidl('\alpha'), xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog,/ylog, ytickformat='logticks_exp'
multiplot
plot, raxis,bcool,xmargin=[8.3,1.7], ystyle=2, xstyle=1 $
      ,charsize=1.5, thick=4,yrange=[1d-2,1d3]*0$
      ,linestyle = 0, ytitle = textoidl('\beta'), xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog,/ylog, ytickformat='logticks_exp'
multiplot
plot, raxis,theta,xmargin=[8.3,1.7], ystyle=2, xstyle=1 $
      ,charsize=1.5, thick=4,yrange=[0.5,2],xtitle='R/AU' $ 
      ,linestyle = 0, ytitle = textoidl('\theta'), xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog
;multiplot
;plot, raxis,1d0/(2d0*!dpi*rate),xmargin=[8.3,1.7], ystyle=2, xstyle=1 $
;      ,charsize=1.5, thick=4 $
;      ,linestyle = 0, ytitle = textoidl('t_{grow}/P_{orb}'), xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog, ytickformat='logticks_exp'
multiplot,/reset
device,/close

xtitle = tex2idl('$R/AU$') + '!X'
ytitle = tex2idl('$t_{grow}/P_{orb}$') + '!X'

if not keyword_set(title) then title = ''

smallh=0.1
tvisc = 1d0/(2d0*!dpi*smallh^2*alpha)

 loadct, 6,/silent
 color_arr = dindgen(2)*256d0/2.

set_plot,'ps'
file = strcompress('result.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, raxis,1d0/(2d0*!dpi*rate),xmargin=[8,2],ymargin=[3.2,2], ystyle=1, xstyle=1 $
      ,charsize=1.5, thick=4, xrange=xrange, title=textoidl(title), xtitle=xtitle, yrange=yrange,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog,/ylog, ytickformat='logticks_exp', color=color_arr(0) 
;oplot, raxis, tvisc, thick=4, linestyle=1
oplot, raxis, 1d0/(2d0*!dpi*rate_special), thick=4, color=color_arr(1)
;oplot, [min(baxis),max(baxis)], [0,0], linestyle=2, thick=2
oplot, [1,1]*raxis(agrid), [1d-8,1d8], thick=2, linestyle=2
xyouts, raxis(agrid)*1.05, 1d2, textoidl('\alpha viscosity>0.1 (t_{c}\Omega<3)'), charsize=1.5

;oplot, [1,1]*raxis(bgrid), [1d-8,1d8], thick=2, linestyle=2
;xyouts, raxis(bgrid)*1.05, 1d4, textoidl('\beta<10'), charsize=1.5


;oplot, [1,1]*raxis[rgrid], [1d-8,1d8], thick=2, linestyle=3
;xyouts,raxis[rgrid]*1.2,3d-2, textoidl('t_{grow}>10 orbits'), charsize=1.5

if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)

    for j=0, n_elements(label)-1 do begin
        ynew = 10.^(alog10(y0) - dy*j)
        oplot, [x0,x1], [ynew,ynew], thick=4, color=color_arr(j)
      ;  print, ynew
      xyouts, x1, ynew, textoidl(label(j)),charsize=1.5
   endfor

endif
device,/close

ytitle = tex2idl('$kH$ for max. growth') + '!X'

set_plot,'ps'
file = strcompress('result_kmax.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, raxis,kmax,xmargin=[8.,2],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
      ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=xtitle, yrange=yrange,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog, color=color_arr(0)


oplot, raxis, kmax_special, thick=4, color=color_arr(1)

oplot, [1,1]*raxis(agrid), [1d-8,1d8], thick=2, linestyle=2
xyouts, raxis(agrid)*1.05, 0.6, textoidl('\alpha viscosity>0.1 (t_{c}\Omega<3)'), charsize=1.5

if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
     for j=0, n_elements(label)-1 do begin
        ynew = y0 - dy*j
        oplot, [x0,x1], [ynew,ynew], thick=4, color=color_arr(j)
      ;  print, ynew
      xyouts, x1, ynew, textoidl(label(j)),charsize=1.5
   endfor



endif

device,/close

ytitle = tex2idl('$\alpha$') + '!X'

set_plot,'ps'
file = strcompress('result_alpha.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, raxis,alpha,xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
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
