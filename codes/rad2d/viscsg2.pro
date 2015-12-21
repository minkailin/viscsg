;assumes constant alpha (cf. constant  nu) and alpha-beta relation holds (i.e. no external heating)
;dispersion relation is for the real growth rate s 
pro viscsg2, bigQ=bigQ, smallq=smallq, gmma=gmma, xrange=xrange, yrange=yrange, kmax=kmax, marg=marg, bmin=bmin, fixalpha=fixalpha 

!p.font = 0 

ii = dcomplex(0d0, 1d0)

if not keyword_set(smallq) then smallq = 1.5
if not keyword_set(gmma) then gmma = 1.4

if keyword_set(marg) then begin ;set Q such that the system is marginally stable in the limit of negligible viscosity 
   bigQ = 1d0/sqrt(2d0*gmma*(2d0-smallq))
endif

kmin = 0d0
if not keyword_set(kmax) then kmax = 1d2
nk = 81920
kaxis = kmin + (kmax-kmin)*dindgen(nk)/(nk-1d0)
rate_fixb = dblarr(nk) 

if not keyword_set(bmin) then bmin = 1d-4
bmax = 1d4

nb = 256
baxis = 10d0^( alog10(bmin) + alog10(bmax/bmin)*dindgen(nb)/(nb-1d0) )

rate_fixb = dblarr(nk) 
rate      = dblarr(nb)

for j=0, nb-1 do begin
bcool = baxis(j) 

;bcool=0.01

if not keyword_set(fixalpha) then begin
   alpha = 1d0
   alpha/= bcool*smallq^2*(gmma-1d0)

   Qext  = 0d0 
endif else begin
   alpha = fixalpha 
   Qext  = 1d0 - alpha*bcool*smallq^2*(gmma-1d0)
endelse 

for i=0, nk-1 do begin
   k = kaxis(i) 


    c5 = bcool^2
  
    c4 = (1d0/3d0)*bcool*(7d0*alpha*bcool*k^2 + 6d0*Qext)

    c3 = (bigQ*(bcool^2*(12d0 + 4d0*alpha^2*k^4 - 6d0*smallq - 6d0*alpha^2*k^2*smallq^2 + $
          3d0*gmma*k^2*(1d0 + 2d0*alpha^2*smallq^2)) + 14d0*alpha*bcool*k^2*Qext + $
          3d0*Qext^2) - 6d0*bcool^2*abs(k))/(3d0*bigQ)

    c2 = (1d0/(3d0*bigQ))*(bigQ*(8d0*alpha^3*bcool^2*(-1d0 + gmma)*k^4*smallq^2 + $
          3d0*bcool*(8d0 + gmma*k^2 - 4d0*smallq)*Qext + $
          2d0*alpha^2*bcool*k^2*(4d0*k^2 + 3d0*(-1d0 + gmma)*smallq^2)*Qext + $
          alpha*k^2*(3d0*bcool^2*(-2d0*(-2d0 + smallq)*smallq + $
          gmma*(k^2 + 2d0*(-1d0 + smallq)*smallq)) + 7d0*Qext^2)) - $ 
          6d0*bcool*(alpha*bcool*k^2 + 2d0*Qext)*Abs(k))

    c1 =  (1d0/(3d0*bigQ))*(bigQ*Qext*(8d0*alpha^3*bcool*(-1d0 + gmma)*k^4*smallq^2 + $
           3d0*alpha*bcool*k^2*(-2d0*(-2d0 + smallq)*smallq + gmma*(k^2 + 2d0*(-1d0 + smallq)*smallq)) + $ 
           4d0*alpha^2*k^4*Qext - 6d0*(-2d0 + smallq)*Qext) - $ 
           6d0*(2d0*alpha^2*bcool^2*(-1d0 + gmma)*k^2*smallq^2 + 2d0*alpha*bcool*k^2*Qext + $
           Qext^2)*abs(k))

   c0 = -((2d0*alpha*k^2*Qext*(2d0*alpha*bcool*(-1d0 + gmma)*smallq^2 + Qext)*abs(k))/bigQ)

    roots = fz_roots([c0, c1, c2, c3, c4, c5])

  
   rate_fixb(i) = max(real_part(roots))

endfor
;stop

rate(j) = max(rate_fixb,grid)

endfor 

xtitle = tex2idl('$\beta$') + '!X'
ytitle = tex2idl('max $Im(\sigma)$') + '!X' + tex2idl('$/\Omega$') + '!X'

set_plot,'ps'
file = strcompress('viscsg2.ps',/remove_all)
device, filename=file $
        ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
plot, baxis,rate,xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
      ,charsize=1.5, thick=4, xrange=xrange, title=title, xtitle=xtitle, yrange=yrange,$
      linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog
;for j=1, 3 do begin
;   oplot, baxis, rate(j,*), thick=4, linestyle=j
;endfor

oplot, [bmin, bmax], [0.,0.], linestyle=2
device,/close
end
