;inviscid limit
;dispersion relation is for growth rate 
pro inviscsg, Qvals=Qvals, smallq=smallq, gmma=gmma, xrange=xrange, yrange=yrange, kmax=kmax, bmin=bmin, $
              legend=legend, label=label, gsoft=gsoft 

!p.font = 0 

ii = dcomplex(0d0, 1d0)

ncases = n_elements(Qvals)

if not keyword_set(smallq) then smallq = 1.5
if not keyword_set(gmma) then gmma = 1.4

bigQ_marg = 1d0/sqrt(2d0*gmma*(2d0-smallq))
print, 'Q for marginal stability is', bigQ_marg

kmin = 0.1
if not keyword_set(kmax) then kmax = 1d2
nk = 8192
kaxis = 10.^(alog10(kmin) + alog10(kmax/kmin)*dindgen(nk)/(nk-1d0))
rate_fixb = dblarr(nk) 

if not keyword_set(bmin) then bmin = 1d-4
bmax = 1d4

nb = 256
baxis = 10d0^( alog10(bmin) + alog10(bmax/bmin)*dindgen(nb)/(nb-1d0) )

rate_fixb = dblarr(nk) 
rate      = dblarr(ncases,nb)
alpha_fixb= dblarr(nk)
alpha     = dblarr(ncases,nb)  
max_stress= dblarr(ncases)

for m=0, ncases-1 do begin

   for j=0, nb-1 do begin
      bcool = baxis(j) 
      
      for i=0, nk-1 do begin
         k = kaxis(i) 
         
         if not keyword_set(gsoft) then begin
              bigQ = Qvals(m)
         endif else  begin
              bigQ = Qvals(m)*(1d0 + k)
         endelse 

         c3 = bcool
         
         c2 = 1d0
         
         c1 = 2d0*k/bigQ - (2d0*(2d0-smallq) + gmma*k^2)
         c1*=-bcool
         
         c0 = -( 2d0*k/bigQ - 2d0*(2d0-smallq) )
         
         roots = fz_roots([c0, c1, c2, c3])
         
         rate_fixb(i) = max(real_part(roots))
         alpha_fixb(i) = rate_fixb(i)/k^2
         
      endfor
      
      rate(m,j) = max(rate_fixb,grid)
      alpha(m,j) = max(alpha_fixb) 
      
   endfor 
   ;max value of s/k^2 in the limit of zero cooling
   max_stress(m) = 3d0^(3d0/2d0)
   max_stress(m)/= 2d0^(7d0/2d0)*(2d0-smallq)^(3d0/2d0)*bigQ^2
endfor

;theoretical values 
rate_theory = (2d0*(2d0 - smallq)/baxis)^(1./3.) ;max rate when the disk is marginally stable in adia limit 


xtitle = tex2idl('$\beta$') + '!X'

loadct, 6
color_arr = dindgen(ncases)*256d0/ncases

ytitle = tex2idl('max$(s/\Omega)$') + '!X' 

set_plot,'ps'
file = strcompress('inviscsg.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, baxis,rate(0,*),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
      ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=xtitle, yrange=yrange,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog, color=color_arr(0)
for j=1, ncases-1 do begin
   oplot, baxis, rate(j,*), thick=4, color=color_arr(j)
endfor
oplot, baxis, rate_theory, thick=4, linestyle=1
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


ytitle= textoidl('max(S/K^2)')

set_plot,'ps'
file = strcompress('inviscsg_alpha.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, baxis,alpha(0,*),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
      ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=xtitle, yrange=yrange,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog, color=color_arr(0)
oplot, [bmin,bmax],[1,1]*max_stress(0), linestyle=1, color=color_arr(0),thick=4
for j=1, ncases-1 do begin
   oplot, baxis, alpha(j,*), thick=4, color=color_arr(j)
   oplot, [bmin,bmax],[1,1]*max_stress(j), linestyle=1, color=color_arr(j),thick=4
endfor
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


