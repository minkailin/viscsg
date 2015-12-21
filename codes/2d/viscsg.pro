pro viscsg, Qvals=Qvals, smallq=smallq, gmma=gmma, xrange=xrange, yrange=yrange, fixalpha = fixalpha, kmax=kmax, bmin=bmin, $
            legend=legend, label=label

!p.font = 0 

ii = dcomplex(0d0, 1d0)

ncases = n_elements(Qvals)

if not keyword_set(smallq) then smallq = 1.5
if not keyword_set(gmma) then gmma = 1.4

bigQ_marg = 1d0/sqrt(2d0*gmma*(2d0-smallq))
print, 'Q for marginal stability is', bigQ_marg

kmin = 0.0
if not keyword_set(kmax) then kmax = 1d4
nk = 8192
kaxis = kmin + (kmax-kmin)*dindgen(nk)/(nk-1d0)
rate_fixb = dblarr(nk) 

if not keyword_set(bmin) then bmin = 1d-4
bmax = 1d4

nb = 256
baxis = 10d0^( alog10(bmin) + alog10(bmax/bmin)*dindgen(nb)/(nb-1d0) )

rate_fixb = dblarr(nk) 
rate      = dblarr(ncases,nb)
kmax      = dblarr(ncases,nb)

for m=0, ncases-1 do begin
   bigQ = Qvals(m)
   
   for j=0, nb-1 do begin
      bcool = baxis(j) 
 
;bcool = 1d-2

      if not keyword_set(fixalpha) then begin
         alpha = 1d0
         alpha/= bcool*smallq^2*(gmma-1d0)
      endif else begin
         alpha = fixalpha
      endelse
      
      for i=0, nk-1 do begin
         k = kaxis(i) 
         
         c4 = bcool
         
         c3 = ii*(1d0 + (7d0/3d0)*alpha*bcool*k^2)
         
         c2 = 2d0*k/bigQ - ( 2d0*(2d0-smallq) + gmma*k^2 + 4d0*alpha^2*k^4/3d0 ) 
         c2*= bcool   
         c2-= 7d0*alpha*k^2/3d0
         
         c1 = 4d0*alpha^2*k^4/3d0 + 2d0*(2d0-smallq) $ 
              + alpha*bcool*k^2*( 3d0*smallq*(2d0-smallq) + gmma*(k^2 +smallq*(3d0*smallq-4d0)) )
         c1*=-ii
         c1+= 2d0*k*ii*(1d0 + alpha*bcool*k^2)/bigQ 
         
         c0 = smallq*(2d0 - alpha*bcool*k^2*(gmma-1d0)*smallq) - 2d0*k/bigQ 
         c0*= alpha*k^2
 
;          c3 =  1d0 + 7d0*k^2/3d0/smallq^2/(gmma-1d0)
;          c2 =  (4d0/3d0)*k^2/(gmma-1d0)/smallq^2 + (7d0/3d0)
;          c2*= alpha*k^2
;          c1 =  (4d0/3d0)*alpha^2*k^4
;          c0 =  2d0*smallq - k^2 - 2d0*k/bigQ
;          c0*= alpha*k^2     
;          rate_fixb(i) = max(real_part(roots))

         roots = fz_roots([c0, c1, c2, c3, c4])
       
;         roots_freq = real_part(roots) 
;         pure_growth = where(roots_freq eq 0d0)

         rate_fixb(i) = max(imaginary(roots))
 
        
      endfor
      
      rate(m,j) = max(rate_fixb,grid)
      kmax(m,j) = kaxis(grid)
   endfor 
endfor

;theoretical rate in the limit of zero viscosity, and disk is
;marginally stable in the limit beta-> infinity (requires Q = sqrt(1/2*gmma*(2-q)))
rate_theory = (2d0*(2d0 - smallq)/baxis)^(1./3.) 

xtitle = tex2idl('$\beta$') + '!X'
ytitle = tex2idl('max$(s/\Omega)$') + '!X' 


loadct, 6
color_arr = dindgen(ncases)*256d0/ncases

set_plot,'ps'
file = strcompress('viscsg.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, baxis,rate(0,*),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
      ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=xtitle, yrange=yrange,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog, color=color_arr(0)
for j=0, ncases-1 do begin
   oplot, baxis, rate(j,*), thick=4, color=color_arr(j)
endfor
;oplot, baxis, rate_theory, thick=4, linestyle=1
if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
   for j=0, ncases-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, color=color_arr(j)
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
endif
device,/close

ytitle = tex2idl('$kH$ for max. growth') + '!X' 

set_plot,'ps'
file = strcompress('viscsg_kmax.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, baxis,kmax(0,*),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
      ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=xtitle, yrange=yrange,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog, color=color_arr(0)
for j=0, ncases-1 do begin
   oplot, baxis, kmax(j,*), thick=4, color=color_arr(j)
endfor
;oplot, baxis, rate_theory, thick=4, linestyle=1
if keyword_set(legend) then begin
   x0=legend(0)
   x1=legend(1)
   y0=legend(2)
   dy=legend(3)
   for j=0, ncases-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4, color=color_arr(j)
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
endif
device,/close

end
