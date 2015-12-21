FUNCTION logticks_exp, axis, index, value
   ; Determine the base-10 exponent
   exponent   = LONG( ALOG10( value ) )
   ; Construct the tickmark string based on the exponent
;   tickmark = '10!E' + STRTRIM( STRING( exponent ), 2 ) + '!N'
    tickmark = textoidl('10^{' + STRTRIM( STRING( exponent ), 2 )+'}')
   ; Return the formatted tickmark string
   RETURN, tickmark
END

function chebyshev_poly, l, zbar

  lsq = l*l
  t = acos(zbar)
  ltimest = l*t
  T_l = cos(ltimest)
  if(abs(zbar) lt 1d0) then begin 
     dT_l = l*sin(ltimest)/sin(t)
     d2T_l= -lsq*cos(ltimest)/sin(t)^2 + l*cos(t)*sin(ltimest)/sin(t)^3
  endif else begin
     dT_l =lsq
     d2t_l =lsq*(lsq-1d0)/3d0
     if(zbar eq -1d0) then begin 
        dT_l = (-1d0)^(l+1d0)*dT_l
        d2T_l= (-1d0)^(l+2d0)*d2T_l
     endif
  endelse

  return, [t_l, dt_l, d2t_l]
end

pro result, xrange=xrange, yrange=yrange, label=label, legend=legend, bslice=bslice, kslice=kslice, vbc=vbc, plot2d=plot2d, $
            title=title
  
  !p.font = 0

  ii          = dcomplex(0d0,1d0)
  one_third   = 1d0/3d0
  two_thirds  = 2d0/3d0
  four_thirds = 4d0/3d0

  if not keyword_set(bslice) then bslice = [1d0]
  if not keyword_set(kslice) then kslice = [1d0]
  if not keyword_set(vbc) then vbc = 'wall'

;parameters
  params = dblarr(11,1)
  openr,1,'params.dat'
  readf,1,params
  close,1
  Q3d     = params(0,0)       ;reference Q 
  smallq  = params(1,0)
  bgmma   = params(2,0)
  gmma    = params(3,0)
  mu      = params(4,0)
  lambda  = params(5,0)
  tirr    = params(6,0)
  fixalpha= params(7,0)
  nmodes  = fix(params(8,0)) ;spectral grid 
  maximize= params(9,0)      ; are we maximizing growth rates over k?
  bulkv    = params(10,0)     ; bulk visc as mult of dyn visc 

;parameter space in k
  nk    = file_lines('kvals.dat')
  kvals = dblarr(nk)
  openr,1,'kvals.dat'
  readf,1,kvals
  close,1

;parameter space in bcool
  nb    = file_lines('bvals.dat')
  bvals = dblarr(nb)
  openr,1,'bvals.dat'
  readf,1,bvals
  close,1

;corresponding Q space  (Q varies with  alpha, and hence beta, if `gviscosity' invoked) 
  Qarr = dblarr(nb)
  openr,1,'Qvals.dat'
  readf,1, Qarr
  close,1

  loadct, 6,/silent
  color_arr = dindgen(nb)*256d0/nb

;basic state 
  nz    = file_lines('basic.dat')
  basic = dblarr(5,nz)
  openr,1,'basic.dat'
  readf,1, basic
  close,1
  zaxis = basic(0,*)
  dens  = basic(1,*)
  dlogd = basic(2,*)
  csq   = basic(3,*)
  dcsq  = basic(4,*)
  zmax  = zaxis(nz-1)
  bigH  = zmax 

;get surface density 
  surf_density = 2d0*int_tabulated(zaxis,dens)
  Q2d = 4d0*Q3d/surf_density
  print,'Q3d, Q2d', Q3d, Q2d

;viscosity  
  visc = dblarr(nb,nz)
  openr,1,'visc.dat'
  readf,1,visc
  close,1

  dvisc= dblarr(nb,nz)
  openr,1,'dvisc.dat'
  readf,1,dvisc
  close,1 

  xtitle = 'z/H'
  ytitle = textoidl('\rho/\rho_0') 
  set_plot,'ps'
  file = strcompress('basic.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, zaxis,dens,xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=yrange,$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval

 ;temporary: analytical density for bgmma=2 (assume omegaz=omega=cs=1)
  ;acoeff = sqrt(1d0/Q3d)
  ;dens_theory = (1d0 + Q3d)*cos(acoeff*zaxis) - Q3d

  ;temporary: gaussian disk (no sg)
  ;dens_theory = exp(-zaxis^2/2d0)

  ;oplot, zaxis, dens_theory, thick=8, linestyle=1
  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, n_elements(label)-1 do begin
        oplot, [x0,x1], [y0,y0]-dy*j, thick=4,linestyle=j 
        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
  endif
  device,/close
  
;eigenvalues
  if(maximize gt 0d0) then begin
   print, 'growth rates maximized over k, only optimal k output'
   nk = 1
  endif 

  ntot = nb*nk
  eigen=dblarr(3,ntot)
  openr,1,'eigenvalues.dat'
  readf,1,eigen
  close,1

  growth = dblarr(nk,nb)
  freq   = dblarr(nk,nb)
  knum   = dblarr(nk,nb)

  for i=0, nb-1 do begin
     ink = i*nk
     growth(*,i) = eigen(0, ink: ink + nk - 1)
     freq(*,i)   = eigen(1, ink: ink + nk - 1) 
     knum(*,i)   = eigen(2, ink: ink + nk - 1)
  endfor

  if(maximize le 0d0) then begin 

  temp = min(abs(bslice[0] - bvals),bgrid)
  print, 'plotting bvals:'
  print, bvals(bgrid), format='(e9.2)'
  
  if not keyword_set(title) then title = '' 
  xtitle = 'kH'
  ytitle = textoidl('Re(s/\Omega)') 
  set_plot,'ps'
  file = strcompress('growth.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, kvals,growth[*,bgrid],xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=yrange,$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval, /xlog, color=color_arr(0),/ylog, title=title 

  for i=1, n_elements(bslice)-1 do begin
     temp = min(abs(bslice[i] - bvals),bgrid)
     print, bvals(bgrid), format='(e9.2)'
     oplot, kvals,growth[*,bgrid],thick=4,color=color_arr(i);,linestyle=i
  endfor

  if keyword_set(plot2d) then begin
  nlines = file_lines('../2d/viscsg_modes.dat')
  array  = dblarr(n_elements(bslice),nlines)
  openr,1,'../2d/viscsg_modes.dat'
  readf,1,array
  close,1
  nlines = file_lines('../2d/viscsg_modes_kaxis.dat')
  kaxis = dblarr(nlines)
  openr,1,'../2d/viscsg_modes_kaxis.dat'
  readf,1,kaxis
  close,1

 
  for i=0, n_elements(bslice)-1 do begin
  oplot, kaxis, array(i,*), thick=4, color=color_arr(i), linestyle=1
  endfor 
  endif 


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
  endif else begin


  loadct, 6,/silent
  color_arr = dindgen(2)*256d0/2.

  xtitle = textoidl('\beta')
  ytitle = textoidl('max(s/\Omega)')   
  if not keyword_set(title) then title = ''
  set_plot,'ps'
  file = strcompress('growth.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, bvals,growth[0,*],xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=1, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=yrange,title=textoidl(title),$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval, /xlog, /ylog, ytickformat='logticks_exp', color=color_arr(0)

;  oplot, bvals, 0.5/bvals^1.5, thick=4, linestyle=2

  if keyword_set(plot2d) then begin ;plot 2d results. need data output file from 2d code 

  for n=0, 1 do begin 
;  if(n eq 0) then file='2dvisc_soft0d64.dat'
;  if(n eq 1) then file='2dvisc_soft2d56.dat'
;  if(n eq 0) then file='2indvisc_soft0d64.dat'
   if(n eq 0) then file='invisc_Hsg0d64.dat'
   if(n eq 1) then file='invisc_Hsg0d52.dat' 
 
  lines = file_lines(file)
  array = dblarr(7,lines)

  openr,1,file

  readf,1,array
  close,1

  baxis = array(0,*)
  rate  = array(2,*)
  kmax  = array(4,*)
 
 
  oplot, baxis , rate, thick=4, psym=2*n+2, symsize=1.5
  endfor 

  endif 
   
  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, n_elements(label)-1 do begin

        ynew = 10.^(alog10(y0) - dy*j)

        if not keyword_set(plot2d) then  begin
         oplot, [x0,x1], [1,1]*ynew, thick=4, linestyle=j
        endif else begin
          if(j eq 0) then begin
         oplot, [x0,x1], [1,1]*ynew, thick=4, linestyle=j
          endif else begin
          oplot, [x1,x1]*0.9, [1,1]*ynew, thick=4, psym=2*j, symsize=1.5
          endelse
           
        endelse


        xyouts, x1, ynew, textoidl(label(j)),charsize=1.5

     endfor
  endif
  device,/close

  ytitle = textoidl('kH for max. growth rate')
  set_plot,'ps'
  file = strcompress('growth_maxk.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, bvals,knum[0,*],xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=yrange,$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval, /xlog, title=textoidl(title)
  if keyword_set(plot2d) then begin 
     for n=0, 1 do begin
   if(n eq 0) then file='invisc_Hsg0d64.dat'
   if(n eq 1) then file='invisc_Hsg0d52.dat'

  lines = file_lines(file)
  array = dblarr(7,lines)

  openr,1,file

  readf,1,array
  close,1

  baxis = array(0,*)
  rate  = array(2,*)
  kmax  = array(4,*)


  oplot, baxis , kmax, thick=4, psym=2*(n+1), symsize=1.5
  endfor
  endif

  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, n_elements(label)-1 do begin    
        if not keyword_set(plot2d) then  begin  
         oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        endif else begin
              ynew = y0-dy*j
            if(j eq 0) then begin
         oplot, [x0,x1], [1,1]*ynew, thick=4, linestyle=j
          endif else begin
          oplot, [x1,x1]*0.9, [1,1]*ynew, thick=4, psym=2*j, symsize=1.5
          endelse

        endelse 

        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
  endif
  
  device,/close
  endelse 

;eigenvectors 
  ntot = nmodes*nk*nb
  coeff= dblarr(12,ntot)
  openr,1,'eigenvectors.dat'
  readf,1, coeff
  close,1
;spectral coeff
  drho = dcomplexarr(nmodes,nk,nb)
  dpres= dcomplexarr(nmodes,nk,nb)
  dpot = dcomplexarr(nmodes,nk,nb)
  vx   = dcomplexarr(nmodes,nk,nb)
  vy   = dcomplexarr(nmodes,nk,nb)
  vz   = dcomplexarr(nmodes,nk,nb)
  for i=0, nb-1 do begin
     inbnk = i*nk*nmodes
     for j=0, nk-1 do begin
        jnmodes = j*nmodes
        
        beg  = inbnk + jnmodes 
        
        drho(*,j,i) = dcomplex(coeff(0, beg : beg + nmodes-1), $
                               coeff(1, beg : beg + nmodes-1) )
        dpres(*,j,i)= dcomplex(coeff(2, beg : beg + nmodes-1), $
                               coeff(3, beg : beg + nmodes-1) )
        dpot(*,j,i) = dcomplex(coeff(4, beg : beg + nmodes-1), $
                               coeff(5, beg : beg + nmodes-1) )
        vx(*,j,i)   = dcomplex(coeff(6, beg : beg + nmodes-1), $
                               coeff(7, beg : beg + nmodes-1) )
        vy(*,j,i) = dcomplex(coeff(8, beg : beg + nmodes-1), $
                             coeff(9, beg : beg + nmodes-1) )
        vz(*,j,i) = dcomplex(coeff(10, beg : beg + nmodes-1), $
                             coeff(11, beg : beg + nmodes-1) )
        
     endfor
  endfor
  
;perturbations in physical space (for a chosen bcool, k)
  density  = dcomplexarr(nz)
  pressure = dcomplexarr(nz)
  potential= dcomplexarr(nz)
  xvel     = dcomplexarr(nz)
  yvel     = dcomplexarr(nz)
  zvel     = dcomplexarr(nz)
  
  ddensity  = dcomplexarr(nz)
  dpressure = dcomplexarr(nz)
  dpotential= dcomplexarr(nz)
  dxvel     = dcomplexarr(nz)
  dyvel     = dcomplexarr(nz)
  dzvel     = dcomplexarr(nz)
  
  d2density  = dcomplexarr(nz)
  d2pressure = dcomplexarr(nz)
  d2potential= dcomplexarr(nz)
  d2xvel     = dcomplexarr(nz)
  d2yvel     = dcomplexarr(nz)
  d2zvel     = dcomplexarr(nz)

 
  temp = min(abs(bvals-bslice[0]),bgrid)
 
  if(maximize le 0d0) then begin  
  temp = min(abs(kvals-kslice),kgrid)
  k    = kvals(kgrid) 
  endif else begin
  kgrid = 0 
  k = knum(0,bgrid) 
  endelse   

  print, 'eigenfunction plot for k, bcool', k, bvals(bgrid)
  print, 'corresponding eigenvalues are ' , growth(kgrid,bgrid), freq(kgrid,bgrid) 

  for i=0, nz-1 do begin ;physical space
     zbar = zaxis[i]/zmax 
 
     for j=1, nmodes  do begin  ;spectral space 
        
                                ;even fields 
        n = 2d0*(j-1d0)
        res = chebyshev_poly(n, zbar)        
        t_l  = res[0]
        dt_l = res[1]/zmax
        d2t_l= res[2]/zmax^2
        

        density(i)  += t_l*drho(j-1, kgrid, bgrid) 
        pressure(i) += t_l*dpres(j-1, kgrid, bgrid)
        potential(i)+= t_l*dpot(j-1, kgrid,bgrid)
        xvel(i)     += t_l*vx(j-1, kgrid, bgrid)
        yvel(i)     += t_l*vy(j-1, kgrid, bgrid)
        
        ddensity(i)  += dt_l*drho(j-1, kgrid, bgrid) 
        dpressure(i) += dt_l*dpres(j-1, kgrid, bgrid)
        dpotential(i)+= dt_l*dpot(j-1, kgrid,bgrid)
        dxvel(i)     += dt_l*vx(j-1, kgrid, bgrid)
        dyvel(i)     += dt_l*vy(j-1, kgrid, bgrid)

        d2density(i)  += d2t_l*drho(j-1, kgrid, bgrid) 
        d2pressure(i) += d2t_l*dpres(j-1, kgrid, bgrid)
        d2potential(i)+= d2t_l*dpot(j-1, kgrid,bgrid)
        d2xvel(i)     += d2t_l*vx(j-1, kgrid, bgrid)
        d2yvel(i)     += d2t_l*vy(j-1, kgrid, bgrid)

                                ;odd field
        m = n + 1d0
        res = chebyshev_poly(m, zbar)
        t_l  = res[0]
        dt_l = res[1]/zmax
        d2t_l= res[2]/zmax^2
        zvel(i)    +=  t_l*vz(j-1, kgrid, bgrid)
        dzvel(i)   += dt_l*vz(j-1, kgrid, bgrid)
        d2zvel(i)  += d2t_l*vz(j-1, kgrid, bgrid)

     endfor
  endfor 
  
  ;error test 
 
  s     = dcomplex(growth(kgrid,bgrid), freq(kgrid,bgrid))

  bcool = bvals(bgrid)
  Q3d   = Qarr(bgrid)
  theta = tirr 

  alpha = visc(bgrid,*)
  alpha_b = bulkv*alpha
  dalpha= dvisc(bgrid,*)
  dlogrhonu = dblarr(nz)
 
  for i=0, nz-1 do begin 
     if(alpha(i) gt 0d0) then begin
        dlogrhonu(i) = dlogd(i) + dalpha(i)/alpha(i)
     endif else dlogrhonu(i) = dlogd(i) 
  endfor 
  dlogrhonu_b = dlogrhonu

  delta_logrhonu = (1d0+mu)*density + lambda*bgmma*pressure
  delta_logrhonu/= csq
 
  ;density eq 
  test1 = s*density/csq + ii*k*xvel + zvel*dlogd + dzvel
  test1/= max(abs(s*density/csq))

  ;energy eq 
  dHvisc = alpha*smallq^2*delta_logrhonu - 2d0*ii*k*alpha*smallq*yvel
  test2 = s*pressure + gmma*ii*k*xvel*csq/bgmma + csq*dlogd*zvel $
          + gmma*csq*dzvel/bgmma + pressure/bcool - (gmma-1d0)*dHvisc 
  test2-= theta*density/bcool/bgmma 
  test2/= max(abs(s*pressure))

  ;xmom eq 
  dFx = alpha*( d2xvel + dlogrhonu*dxvel - four_thirds*k^2*xvel )
  dFx+= ii*k*alpha*(one_third*dzvel + dlogrhonu*zvel)
  dFx-= alpha_b*k^2*xvel
  dFx+= ii*k*alpha_b*dzvel 
  test3 = s*xvel + ii*k*pressure + ii*k*potential - 2d0*yvel - dFx
  test3/= max(abs(pressure))/bigH 

  ;ymom eq 
  dFy = -ii*k*alpha*smallq*delta_logrhonu
  dFy+= alpha*( d2yvel + dlogrhonu*dyvel - k^2*yvel)
  test4 = s*yvel + (2d0-smallq)*xvel - dFy
  test4/= max(abs(pressure))/bigH

  ;zmom eq 
  dFz = ii*k*alpha*( one_third*dxvel - two_thirds*dlogrhonu*xvel)
  dFz+= alpha*( four_thirds*d2zvel + four_thirds*dlogrhonu*dzvel - k^2*zvel)
  dFz+= alpha_b*(d2zvel + dlogrhonu_b*dzvel)
  dFz+= ii*alpha_b*k*(dxvel + dlogrhonu_b*xvel)
  test5 = s*zvel + dpressure + dlogd*pressure - dlogd*density + dpotential - dFz
  test5/= max(abs(dpressure)) 

  ;poisson eq
  test6 = d2potential - k*k*potential - dens*density/csq/Q3d
  test6/= max(abs(dens*density/csq/Q3d))

  test1 = max(abs(test1(0:nz-2)));/max(abs(s*density))
  test2 = max(abs(test2(0:nz-2)));/max(abs(s*pressure))
  test3 = max(abs(test3(0:nz-2)));/max(abs(s*xvel))
  test4 = max(abs(test4(0:nz-2)));/max(abs(s*yvel))
  test5 = max(abs(test5(0:nz-2)));/max(abs(s*zvel))
  test6 = max(abs(test6(0:nz-2)));/max(abs(k*k*potential))

  internal_error = max([test1,test2,test3,test4,test5,test6])
  print, 'internal err', internal_error

  ;boundary error
  if(vbc eq 'wall') then begin
  bound1 = dpotential(nz-1) + k*potential(nz-1)  
  bound1/= max(abs(potential))/bigH 
  endif
  if(vbc eq 'free') then begin
  bound1 =s*(dpotential(nz-1) + k*potential(nz-1))  + dens(nz-1)*zvel(nz-1)/Q3d 
  bound1/= max(abs(s*potential))/bigH 
  endif

  bound2 = dxvel(nz-1)
  bound2/= max(abs(xvel))/bigH + 1d-16
  bound3 = dyvel(nz-1)
  bound3/= max(abs(yvel))/bigH + 1d-16

  if(vbc eq 'wall') then begin
  bound4 = zvel(nz-1)
  bound4/= max(abs(zvel)) + 1d-16 
  endif 
  if(vbc eq 'free') then begin
  bound4 = s*pressure(nz-1) + csq(nz-1)*dlogd(nz-1)*zvel(nz-1)
  bound4/= max(abs(s*pressure))
  endif 

  boundary_error = max(abs([bound1,bound2,bound3,bound4]))
  print, 'boundary error', boundary_error

  ;integral test for inviscid channel mode in vert iso disk 
  if( (max(alpha) eq 0d0) and (Bgmma eq 1d0) and (k eq 0d0) and (vbc eq 'wall')) then begin 
  inttest = s*s*int_tabulated(zaxis,dens*abs(zvel)^2)
  bot = inttest 
  inttest-= (1d0/Q3d)*int_tabulated(zaxis,dens*dens*abs(zvel)^2)
  temp = dcomplex(deriv(zaxis,real_part(dens*zvel)), deriv(zaxis,imaginary(dens*zvel)))
  inttest+= int_tabulated(zaxis, abs(temp)^2/dens )
  inttest+= (bcool*s*gmma/(1d0+bcool*s) - 1d0)*int_tabulated(zaxis,dens*abs(dzvel)^2)
  temp = dlogd*dens*zvel*conj(dzvel) 
  inttest+= (bcool*s/(1d0+bcool*s) - 1d0)*dcomplex(int_tabulated(zaxis,real_part(temp)),int_tabulated(zaxis,imaginary(temp)))

  print,'integral test for inviscid channel mode', abs(inttest)/abs(bot) 
  endif 





  density *= dens/csq
  density *= conj(density(0)) 
  density /= max(abs(density)) 
  xtitle = 'z/H'
  ytitle = textoidl('\delta\rho')
  set_plot,'ps'
  file = strcompress('eigenvec_dens.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, zaxis,real_part(density),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=[-1,1],$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval
  oplot, zaxis, imaginary(density), thick=4, linestyle=1 
  device,/close

  pressure *= dens 
  pressure *= conj(pressure(0))
  pressure /= max(abs(pressure))
  ytitle = textoidl('\deltaP')
  set_plot,'ps'
  file = strcompress('eigenvec_pres.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, zaxis,real_part(pressure),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=[-1,1],$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval
  oplot, zaxis, imaginary(pressure), thick=4, linestyle=1
  device,/close

  zvel *= conj(zvel(0))
  zvel /= max(abs(zvel))
  ytitle = textoidl('\deltav_z')
  set_plot,'ps'
  file = strcompress('eigenvec_vz.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, zaxis,real_part(zvel),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=[-1,1],$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval
  oplot, zaxis, imaginary(zvel), thick=4, linestyle=1
  device,/close

end 

