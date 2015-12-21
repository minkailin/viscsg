;dispersion relation for the real growth rate s
;for viscosity prescription nu propto Sigma^mu Pressure^lambda
pro viscsg_modes, bigQ=bigQ, smallq=smallq, gmma=gmma, mu=mu, lambda=lambda, xrange=xrange, yrange=yrange, fixalpha = fixalpha, kmin=kmin, kmax=kmax, $
            legend=legend, label=label, nk=nk, gsoft=gsoft, bvals=bvals, title=title, tirr = tirr, bulkv = bulkv, maximize=maximize, nb=nb, gvisc=gvisc 

!p.font = 0 

if not keyword_set(smallq) then smallq = 1.5
if not keyword_set(gmma) then gmma = 1.4
if not keyword_set(mu) then mu = 0d0
if not keyword_set(lambda) then lambda = 0d0 
if not keyword_set(tirr) then tirr = 0d0 
if not keyword_set(bulkv) then bulkv =0d0 

if not keyword_set(gsoft) then begin
   grav3d = 0d0
   bigQ_marg = 1d0/sqrt(2d0*gmma*(2d0-smallq))
endif else begin
   grav3d = gsoft 
   kcrit = -(1./6.) + gmma/( $
           3.*2.^(2./3.)*(432.*gmma^2. - 2.*gmma^3. - 216.*gmma^2.*smallq + $
                          sqrt(-4.*gmma^6. + (432.*gmma^2. - 2.*gmma^3. - 216.*gmma^2*smallq)^2.))^( $
           1./3.)) + (432.*gmma^2. - 2.*gmma^3. - 216.*gmma^2*smallq + $
                      sqrt(-4.*gmma^6. + (432.*gmma^2. - 2.*gmma^3. - 216.*gmma^2*smallq)^2.))^( $
           1./3.)/(6.*2.^(1./3.)*gmma)
   bigQ_marg = 1d0/gmma/kcrit/(1d0+kcrit)^2.
endelse
bigQ_fix = bigQ*bigQ_marg
print, 'bigQ=', bigQ_fix

if not keyword_set(kmin) then kmin = 1d-2
if not keyword_set(kmax) then kmax = 1d2
if not keyword_set(nk) then nk =8192
dlogk = alog10(kmax/kmin)/(nk-1d0)
kaxis = 10.^(alog10(kmin) + dindgen(nk)*dlogk)

;if keyword_set(gvisc) then begin
;   if not keyword_set(fixalpha) then  begin ;re-set bmin because we need alpha < 1 
;      bvals[0] = 1d0 - tirr
;      bvals[0] =  bvals[0]/( (gmma-1d0)*smallq^2d0 )
;      print, 'reset bmin to', bvals[0]
;   endif
;endif

if keyword_set(maximize) then begin
   if not keyword_set(nb) then nb   = 64
   bmin = bvals[0]
   bmax = bvals[1]
   bvals = 10.^(alog10(bmin)+alog10(bmax/bmin)*dindgen(nb)/(nb-1d0))
endif
nb = n_elements(bvals)

roots = dblarr(nb, 4, nk) 
roots_output = dblarr(nb,nk) 
maxroot_fixb = dblarr(nk)
maxroot = dblarr(nb) 
maxk    = dblarr(nb)

freqs = dblarr(nb, 4, nk)

roots_asymp = dblarr(nb)
roots_gammie_largek = dblarr(nb, nk) 
roots_gammie_smallk = dblarr(nb, nk)

for j=0, nb-1 do begin
   bcool = bvals(j)
   
   if not keyword_set(fixalpha) then begin
         theta = tirr
         alpha = 1d0 - theta
         alpha/= (gmma-1d0)*bcool*smallq^2
   endif else begin             ; alpha, theta bcool are independent 
      alpha = fixalpha
      theta = tirr 
   endelse
   bvisc = bulkv*alpha 

   if keyword_set(gvisc) then begin ;set Q according to alpha 
      bigQ_fix = bigQ_marg/sqrt(alpha)
   endif

   eta = 1d0 - alpha*bcool*smallq^2*(gmma-1d0)*lambda 
   
   ;asymptotic rate for lambda=0, large k 
   
   a = bcool*(4d0*alpha/3d0 + bvisc)
   b = 4d0*alpha/3d0 + bvisc + gmma*bcool
   c = theta - alpha*smallq^2*(1d0+mu)*(gmma-1d0)*bcool
   roots_asymp(j) = -b + sqrt(b^2 - 4d0*a*c)
   roots_asymp(j)/= 2d0*a 
   
   for i=0, nk-1 do begin
      k = kaxis(i) 
      
      bigQ = bigQ_fix*(1d0 + abs(k)*grav3d)
      
   ;aspymtotic rate for mu=-1, lambda=0, small and large 
      roots_gammie_smallk(j,i) = alpha*k^2/(2d0*(2d0-smallq))
      roots_gammie_smallk(j,i)*= 2d0*k/bigQ  - theta*k^2 

      roots_gammie_largek(j,i) = 2d0/bigQ/k - theta 
      roots_gammie_largek(j,i)/= 4d0*alpha/3d0 + bvisc + gmma*bcool 
      

      c0 = (1./bigQ)*alpha*k^2.*(2.*eta*(-k + bigQ*(1. + mu)*smallq) - $
                                 alpha*bcool*(-1. + gmma)*smallq^2*(4.*k*lambda + $
                                                                    bigQ*(1. + mu)*(k^2. - 2.*lambda*smallq)) + $
                                 bigQ*(k^2. + 2.*lambda*smallq)*theta)

      c1 = (1./(3.*bigQ))*(-6.*k*(eta + alpha*bcool*k^2) + $
                           bigQ*(eta*(12. + 4.*alpha^2.*k^4. + 3.*alpha*bvisc*k^4. - 6.*smallq) + $
                                 k^2.*(8.*alpha^3.*bcool*(-1. + gmma)*k^2.*lambda*smallq^2. + $
                                       6.*alpha^2*bcool*bvisc*(-1. + gmma)*k^2.*lambda*smallq^2. + $
                                       3.*alpha*bcool*(-(3. + mu)*(-2. + smallq)*smallq + $
                                                       gmma*(k^2. + smallq*(-4. + 2.*lambda + (3. + mu)*smallq))) + $
                                       3.*theta)))

      c2 = (1./(3.*bigQ))*(bigQ*(7.*alpha + 3.*bvisc)*eta*k^2 + $
                           bcool*(-6.*k + $
                                  bigQ*(12. + 4.*alpha^2.*k^4. + 3.*alpha*bvisc*k^4. - 6.*smallq - $
                                        6.*alpha^2.*k^2.*lambda*smallq^2. + $
                                        3.*gmma*k^2.*(1. + 2.*alpha^2.*lambda*smallq^2.))))

      c3 = eta + (1./3.)*bcool*(7.*alpha + 3.*bvisc)*k^2
      
      c4 = bcool 
      
      solve = fz_roots([c0,c1,c2,c3,c4])
      
      roots(j,*,i) = real_part(solve)
      freqs(j,*,i) = imaginary(solve)
      
      maxroot_fixb(i)  = max(real_part(solve))
      roots_output(j,i)= maxroot_fixb(i)
      
   endfor
   maxroot(j) = max(maxroot_fixb,kgrid)
   maxk(j)    = kaxis(kgrid) 
endfor

;temp = min(abs(maxk-1d0),bgrid)
;print, 'beta for max k = 1 is', bvals(bgrid)


openw,1,'viscsg_modes.dat'
printf,1,roots_output
close,1
openw,1,'viscsg_modes_kaxis.dat'
for i=0, nk-1 do begin
   printf,1,kaxis(i)
endfor 
close,1

 if not keyword_set(title) then title = ''

  loadct, 6,/silent
  color_arr = dindgen(nb)*256d0/nb

  if not keyword_set(maximize) then begin
     
     if not keyword_set(yrange) then begin
        yrange1=[0d0,max(roots)]
     endif else begin
        yrange1=yrange
     endelse
     
     xtitle = textoidl('kH') + '!X'
     ytitle = textoidl('Re(s/\Omega)') + '!X'
;xtitle = 'Dimensionless wavenumber '+textoidl('kH')
;ytitle = 'Growth rate in '+textoidl('\Omega')

     set_plot,'ps'
     file = strcompress('viscsg_modes.ps',/remove_all)
     device, filename=file $
             ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
     plot, kaxis, roots(0,0,*),xmargin=[8.,2],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
           ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=yrange1,$
           linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2, color=color_arr(0), /xlog , /ylog, $
     title = textoidl(title)
     oplot, [kmin, kmax], roots_asymp(0)*[1,1], thick=1, linestyle=2
;     oplot, kaxis, roots_gammie_largek(0,*), thick=4,color=color_arr(0), linestyle=1
;     oplot, kaxis, roots_gammie_smallk(0,*), thick=4,color=color_arr(0), linestyle=2
     for j=1, 3 do begin
        oplot, kaxis, roots(0,j,*), thick=4,color=color_arr(0)
     endfor
     
     for i=1,nb-1 do begin
        oplot, [kmin, kmax], roots_asymp(i)*[1,1], thick=1, linestyle=2, color=color_arr(i)
        for j=0, 3 do begin
           oplot, kaxis, roots(i,j,*), thick=4,color=color_arr(i)
        endfor
;        oplot, kaxis, roots_gammie_largek(i,*), thick=4,color=color_arr(i), linestyle=1
;        oplot, kaxis, roots_gammie_smallk(i,*), thick=4,color=color_arr(i), linestyle=2
     endfor
     
     if keyword_set(legend) then begin
        x0=legend(0)
        x1=legend(1)
        y0=legend(2)
        dy=legend(3)
        for j=0, n_elements(label)-1 do begin
           ynew = 10.^(alog10(y0) - dy*j)
           oplot, [x0,x1], [ynew,ynew], thick=4, color=color_arr(j)
           xyouts, x1, ynew, textoidl(label(j)),charsize=1.5
        endfor
     endif
     device,/close
     
  endif else begin
     
     xtitle=textoidl('\beta')
     ytitle=textoidl('max(s/\Omega)')
     set_plot,'ps'
     file = strcompress('viscsg_modes_max.ps',/remove_all)
     device, filename=file $
             ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
     plot, bvals, maxroot,xmargin=[8.,2],ymargin=[3.2,1.8], ystyle=1, xstyle=1 $
           ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=yrange,$
           linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2, color=color_arr(0), /xlog ;, /ylog ;, $
     device,/close
     
     ytitle=textoidl('kH for max. growth')
     set_plot,'ps'
     file = strcompress('viscsg_modes_kmax.ps',/remove_all)
     device, filename=file $
             ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
     plot, bvals, maxk,xmargin=[8.,2],ymargin=[3.2,1.8], ystyle=1, xstyle=1 $
           ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=[kmin,kmax],$
           linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2, color=color_arr(0), /xlog ;, $
     device,/close
  endelse
  
   if not keyword_set(maximize) then begin
     xtitle = textoidl('kH') + '!X'
     ytitle = textoidl('Im(s/\Omega)') + '!X'

;if not keyword_set(yrange) then begin
;   yrange1=[min(freqs(where(roots ge 0d0))),max(freqs(where(roots ge 0d0)))]
;endif else begin
;   yrange1=yrange
;endelse

     set_plot,'ps'
     file = strcompress('viscsg_modes_freq.ps',/remove_all)
     device, filename=file $
             ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
     plot, kaxis, freqs(0,0,*),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
           ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=yrange,$
           linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2, color=color_arr(0), psym=2, symsize=0.3,/xlog ;, $
                                ;/nodata;, title=textoidl(title)
     for j=1, 3 do begin
        growth = roots(0,j,*)
        oscill = freqs(0,j,*)
        filter = where(growth ge 0d0)

        if(min(filter) ge 0) then oplot, kaxis(filter), oscill(filter), thick=4,color=color_arr(0),psym=2, symsize=0.3
     endfor

     for i=0,nb-1 do begin
        for j=0, 3 do begin
           growth = roots(i,j,*)
           oscill = freqs(i,j,*)
           filter = where(growth ge 0d0)

           if(min(filter) ge 0) then oplot, kaxis(filter), oscill(filter), thick=4,color=color_arr(i),psym=2, symsize=0.3
        endfor
     endfor

       if keyword_set(legend) then begin
        x0=legend(0)
        x1=legend(1)
        y0=legend(2)
        dy=legend(3)
        for j=0, n_elements(label)-1 do begin
           oplot, [x0,x1], [y0,y0]-dy*j, thick=4, color=color_arr(j)
           xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
        endfor
     endif
     device,/close
  endif

end
