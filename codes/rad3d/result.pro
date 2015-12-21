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

FUNCTION logticks_exp, axis, index, value
   exponent   = LONG( ALOG10( value ) )
    tickmark = textoidl('10^{' + STRTRIM( STRING( exponent ), 2 )+'}')
   RETURN, tickmark
END

pro result, xrange=xrange, yrange=yrange, label=label, legend=legend, rslice=rslice, kslice=kslice, vbc=vbc, plot2d=plot2d, $
            title=title  
  
  !p.font = 0

  bigG = 6.67d-8 
  mstar = 1.99d33
  au=1.50d13 
  stefan = 5.67d-5
  mean_mu = 2.3
  kboltz = 1.38d-16
  mh=1.66d-24
  rstar = kboltz/mh 

  ii          = dcomplex(0d0,1d0)
  one_third   = 1d0/3d0
  two_thirds  = 2d0/3d0
  four_thirds = 4d0/3d0

  if not keyword_set(rslice) then rslice = [100d0]
  if not keyword_set(kslice) then kslice = [1d0]
  if not keyword_set(vbc) then vbc = 'wall'

;parameters
  params = dblarr(10,1)
  openr,1,'params.dat'
  readf,1,params
  close,1
  
  smallq  = params(0,0)
  gmma    = params(1,0)
  mu      = params(2,0)
  lambda  = params(3,0)
  tirr    = params(4,0)
  nmodes  = fix(params(5,0)) ;spectral grid size
  maximize= params(6,0)      ;are we maximizing growth rates over k?
  bulkv   = params(7,0)      ;bulk visc as mult of dyn visc 
  thincool= params(8,0)      ; rad diff or opt thin cool?
  smallb  = params(9,0)      ; opacity law 

  if (smallb eq 0d0) then opacity0 = 0.24
  if (smallb eq 2d0) then opacity0 = 5d-4 


  ;2d basic state 
  nr      = file_lines('basic2d.dat')
  basic2d = dblarr(9, nr) 
  openr,1,'basic2d.dat'
  readf,1,basic2d
  close,1
  
  raxis    = basic2d(0,*)
  Qarr     = basic2d(3,*)
  alpha_arr= basic2d(4,*)
  bcool2d  = basic2d(5,*)
  theta2d  = basic2d(6,*)
  temp2d   = basic2d(7,*)
  tau2d    = basic2d(8,*)

  res = min(abs(alpha_arr - 0.1),agrid)
  print, 'bcool at alpha=0.1', bcool2d(agrid)

  if(nr gt 1) then begin
     set_plot,'ps'
     file = strcompress('result_basic2d.ps',/remove_all)
     device, filename=file $
             ,bits_per_pixel=8,xsize=8, ysize=6,xoffset=0,yoffset=0,/inches,/color
     multiplot,/reset
     multiplot,[1,4],margins=[0.12,0.1,0.05,0.05],rowspacing=0.1
     plot, raxis,Qarr,ymargin=[5,0], ystyle=1, xstyle=1 $
           ,charsize=1.5, thick=4,yrange=[1,10],$
           linestyle = 0, ytitle = 'Q', xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog
     multiplot
     plot, raxis,alpha_arr, ystyle=2, xstyle=1 $
           ,charsize=1.5, thick=4 $
           ,linestyle = 0, ytitle = textoidl('\alpha'), xtickinterval=xtickinterval, ytickinterval=ytickinterval,$
           charthick=2,/xlog,/ylog, ytickformat='logticks_exp';, xtitle='R/AU'
     multiplot
     plot, raxis, bcool2d, ystyle=2, xstyle=1 $
           ,charsize=1.5, thick=4 $
           ,linestyle = 0, ytitle = textoidl('\beta'), xtickinterval=xtickinterval, ytickinterval=ytickinterval,$
           charthick=2,/xlog,/ylog, ytickformat='logticks_exp'
     multiplot
     plot, raxis, tau2d,ymargin=[5,0], ystyle=2, xstyle=1 $
           ,charsize=1.5, thick=4, xtitle='R/AU',$
           linestyle = 0, ytitle = textoidl('\tau'), xtickinterval=xtickinterval, ytickinterval=ytickinterval,charthick=2,/xlog,/ylog 
     multiplot,/reset
     device,/close
  endif

  loadct, 6,/silent
  nsample    = n_elements(rslice) 
  color_arr = dindgen(nsample)*256d0/nsample 

; 3d basic state. each radius has a different veritcal structure 
  nztot   = file_lines('basic3d.dat') ;nztot = nz*nr 
  nz      = nztot/nr 

  basic3d = dblarr(11,nztot)
  openr,1,'basic3d.dat'
  readf,1, basic3d
  close,1
  
  res  = min(abs(raxis - rslice(0)),rgrid)
  rbeg = (rgrid)*nz
  rend = rbeg + nz - 1

  zaxis = basic3d(0,rbeg:rend)
  dens  = basic3d(1,rbeg:rend)
  temp  = basic3d(5,rbeg:rend) 

  xtitle = 'z/H(R)'
  ytitle = textoidl('T/T_0') 
  set_plot,'ps'
  file = strcompress('result_basic3d.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, zaxis,temp,xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=yrange,$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval, color = color_arr(0) 
  for i=1, nsample-1 do begin
  temp = min(abs(raxis - rslice(i)),rgrid)
  rbeg = (rgrid)*nz
  rend = rbeg + nz - 1

  zaxis = basic3d(0,rbeg:rend)
  dens  = basic3d(1,rbeg:rend)
  temp  = basic3d(5,rbeg:rend)     

  oplot, zaxis, temp, thick=4, color=color_arr(i) 
  endfor 

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


;optical depth profile 

  res    = min(abs(raxis - rslice(0)),rgrid)
  radius = raxis(rgrid)
  Q3d    = basic2d(2,rgrid)
  T0     = basic2d(7,rgrid)
  omega  = sqrt(bigG*mstar/(radius*au)^3)
  rhomid = omega^2/(4d0*!dpi*bigG*Q3d)
  cs0    = sqrt(rstar*T0/mean_mu)

  rbeg = (rgrid)*nz
  rend = rbeg + nz - 1

  zaxis = basic3d(0,rbeg:rend)
  dens  = basic3d(1,rbeg:rend)
  temp  = basic3d(5,rbeg:rend) 
  
  tau0 = opacity0*T0^2*(cs0/omega)*rhomid
 
  tauz = dblarr(nz)

  for j=0, nz-3 do begin
     tauz(j) = int_tabulated(zaxis(j:nz-1), temp(j:nz-1)^2*dens(j:nz-1))
  endfor
  tauz *= tau0 

  ytitle = textoidl('\tau') 
  set_plot,'ps'
  file = strcompress('result_basic3d_tau.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, zaxis(0:nz-3),tauz(0:nz-3),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=[1d-1,1d5],$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval, color = color_arr(0),/ylog 
  for i=1, nsample-1 do begin
     res    = min(abs(raxis - rslice(i)),rgrid)
     radius = raxis(rgrid)
     Q3d    = basic2d(2,rgrid)
     T0     = basic2d(7,rgrid)
     omega  = sqrt(bigG*mstar/(radius*au)^3)
     rhomid = omega^2/(4d0*!dpi*bigG*Q3d)
     cs0    = sqrt(rstar*T0/mean_mu)
     
     rbeg = (rgrid)*nz
     rend = rbeg + nz - 1
     
     zaxis = basic3d(0,rbeg:rend)
     dens  = basic3d(1,rbeg:rend)
     temp  = basic3d(5,rbeg:rend) 
     
     tau0 = opacity0*T0^2*(cs0/omega)*rhomid
    
     for j=0, nz-3 do begin
        tauz(j) = int_tabulated(zaxis(j:nz-1), temp(j:nz-1)^2*dens(j:nz-1))
     endfor
     tauz *= tau0 
     
     oplot, zaxis(0:nz-3), tauz(0:nz-3), thick=4, color=color_arr(i) 
  endfor 

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



;parameter space in k
  nk    = file_lines('kvals.dat')
  kvals = dblarr(nk)
  openr,1,'kvals.dat'
  readf,1,kvals
  close,1

;eigenvalues
  if(maximize gt 0d0) then begin
   print, 'growth rates maximized over k, only optimal k output'
   nk = 1
  endif 

  ntot = nr*nk
  eigen=dblarr(3,ntot)
  openr,1,'eigenvalues.dat'
  readf,1,eigen
  close,1

  growth = dblarr(nk,nr)
  freq   = dblarr(nk,nr)
  knum   = dblarr(nk,nr)

  for i=0, nr-1 do begin
     ink = i*nk
     growth(*,i) = eigen(0, ink: ink + nk - 1)
     freq(*,i)   = eigen(1, ink: ink + nk - 1) 
     knum(*,i)   = eigen(2, ink: ink + nk - 1)
  endfor

  if(maximize le 0d0) then begin 

  temp = min(abs(rslice[0] - raxis),rgrid)
  print, 'plotting radii:'
  print, raxis(rgrid), format='(e9.2)'
  
  if not keyword_set(title) then title = '' 
  xtitle = 'kH'
  ytitle = textoidl('Re(s/\Omega)') 
  set_plot,'ps'
  file = strcompress('result_growth.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, kvals,growth[*,rgrid],xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=yrange,$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval, color=color_arr(0),/ylog, title=title, ytickformat='logticks_exp' 

  for i=1, nsample-1 do begin
     temp = min(abs(rslice[i] - raxis),rgrid)
     print, raxis(rgrid), format='(e9.2)'
     oplot, kvals,growth[*,rgrid],thick=4,color=color_arr(i)
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

  xtitle = textoidl('R/AU')
  ytitle = textoidl('t_{grow}/P_{orb}')   
  if not keyword_set(title) then title = ''
  set_plot,'ps'
  file = strcompress('result_growth.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, raxis,1d0/(2d0*!dpi*growth[0,*]),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=1, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=yrange,title=textoidl(title),$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval, /xlog, /ylog, color=color_arr(0), ytickformat='logticks_exp'
  
  oplot, [1,1]*raxis(agrid), [1d-8,1d8], thick=2, linestyle=2
  xyouts, raxis(agrid)*1.05, 1d4, textoidl('\alpha>0.1'), charsize=1.5
   xyouts, raxis(agrid)*1.05, 4d3, textoidl('(t_c\Omega<4)'), charsize=1.5

  if keyword_set(plot2d) then begin ;plot 2d results. need data output file from 2d code 

  lines = file_lines('../rad2d/output.dat')
  array = dblarr(11,lines)

  openr,1,'../rad2d/output.dat'
  readf,1,array
  close,1

  raxis2d = array(0,*)
  rate  = array(2,*)
  kmax  = array(4,*)
 
  oplot, raxis2d , 1d0/(2d0*!dpi*rate), thick=4, color=color_arr(1)

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
;            if (j eq n_elements(label)-1) then begin
;               oplot, [x1,x1]*0.95, [1,1]*ynew, thick=4, psym=2, symsize=1.5
;         endif else begin
         oplot, [x0,x1], [1,1]*ynew, thick=4, color=color_arr(j)
;         endelse
        endelse


        xyouts, x1, ynew, textoidl(label(j)),charsize=1.5

     endfor
  endif
  device,/close

  ytitle = textoidl('kH for max. growth rate')
  set_plot,'ps'
  file = strcompress('result_growth_maxk.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, raxis,knum[0,*],xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=yrange,$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval, /xlog, title=textoidl(title), color=color_arr(0)
  if keyword_set(plot2d) then oplot, raxis2d, kmax, thick=4, color=color_arr(1)

  oplot, [1,1]*raxis(agrid), [1d-8,1d8], thick=2, linestyle=2
  xyouts, raxis(agrid)*1.05, 0.4, textoidl('\alpha>0.1'), charsize=1.5
   xyouts, raxis(agrid)*1.05, 0.3, textoidl('(t_c\Omega<4)'), charsize=1.5



  if keyword_set(legend) then begin
     x0=legend(0)
     x1=legend(1)
     y0=legend(2)
     dy=legend(3)
     for j=0, n_elements(label)-1 do begin    
        if not keyword_set(plot2d) then  begin  
         oplot, [x0,x1], [y0,y0]-dy*j, thick=4, linestyle=j
        endif else begin
;         if (j eq n_elements(label)-1) then begin
;         oplot, [x1,x1]*0.9, [y0,y0]-dy*j, thick=4, psym=2, symsize=1.5
;         endif else begin
         oplot, [x0,x1], [y0,y0]-dy*j, thick=4, color=color_arr(j)
;         endelse
        endelse 

        xyouts, x1, y0-dy*j,textoidl(label(j)),charsize=1.5
     endfor
  endif
  
  device,/close
  endelse 



;eigenvectors 
  ntot = nmodes*nk*nr
  coeff= dblarr(12,ntot)
  openr,1,'eigenvectors.dat'
  readf,1, coeff
  close,1
;spectral coeff
  dpres= dcomplexarr(nmodes,nk,nr)
  dtmp = dcomplexarr(nmodes,nk,nr)
  dpot = dcomplexarr(nmodes,nk,nr)
  vx   = dcomplexarr(nmodes,nk,nr)
  vy   = dcomplexarr(nmodes,nk,nr)
  vz   = dcomplexarr(nmodes,nk,nr)
  for i=0, nr-1 do begin
     inbnk = i*nk*nmodes
     for j=0, nk-1 do begin
        jnmodes = j*nmodes
        
        beg  = inbnk + jnmodes 
        
        dpres(*,j,i) = dcomplex(coeff(0, beg : beg + nmodes-1), $
                               coeff(1, beg : beg + nmodes-1) )
        dtmp(*,j,i)= dcomplex(coeff(2, beg : beg + nmodes-1), $
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
  temp     = dcomplexarr(nz) 
  potential= dcomplexarr(nz)
  xvel     = dcomplexarr(nz)
  yvel     = dcomplexarr(nz)
  zvel     = dcomplexarr(nz)
  
  ddensity  = dcomplexarr(nz)
  dpressure = dcomplexarr(nz)
  dtemp     = dcomplexarr(nz)
  dpotential= dcomplexarr(nz)
  dxvel     = dcomplexarr(nz)
  dyvel     = dcomplexarr(nz)
  dzvel     = dcomplexarr(nz)
  
  d2density  = dcomplexarr(nz)
  d2pressure = dcomplexarr(nz)
  d2temp     = dcomplexarr(nz) 
  d2potential= dcomplexarr(nz)
  d2xvel     = dcomplexarr(nz)
  d2yvel     = dcomplexarr(nz)
  d2zvel     = dcomplexarr(nz)

  res = min(abs(raxis-rslice[0]),rgrid)
 
  if(maximize le 0d0) then begin  
  res = min(abs(kvals-kslice[0]),kgrid)
  k    = kvals(kgrid) 
  endif else begin
  kgrid = 0 
  k = knum(0,rgrid) 
  endelse   

  print, 'eigenfunction plot for k, radius', k, raxis(rgrid)
  print, 'corresponding eigenvalues are ' , growth(kgrid,rgrid), freq(kgrid,rgrid) 

  rbeg = (rgrid)*nz
  rend = rbeg + nz - 1

  zaxis  = basic3d(0,rbeg:rend)
  d      = basic3d(1,rbeg:rend)
  dlogd  = basic3d(2,rbeg:rend) 
  p      = basic3d(3,rbeg:rend)
  dlogp  = basic3d(4,rbeg:rend)
  t      = basic3d(5,rbeg:rend)
  dlogt  = basic3d(6,rbeg:rend)
  d2logt = basic3d(7,rbeg:rend)
  sg     = basic3d(8,rbeg:rend)
  dsg    = basic3d(9,rbeg:rend)
  fz     = basic3d(10,rbeg:rend)

  zmax  = basic2d(1,rgrid) 

  for i=0, nz-1 do begin ;physical space
     zbar = zaxis[i]/zmax 
 
     for j=1, nmodes  do begin  ;spectral space 
        
                                ;even fields 
        n = 2d0*(j-1d0)
        res = chebyshev_poly(n, zbar)        
        t_l  = res[0]
        dt_l = res[1]/zmax
        d2t_l= res[2]/zmax^2
        

        pressure(i) += t_l*dpres(j-1, kgrid, rgrid)
        temp(i)     += t_l*dtmp(j-1, kgrid, rgrid)
        potential(i)+= t_l*dpot(j-1, kgrid,rgrid)
        xvel(i)     += t_l*vx(j-1, kgrid, rgrid)
        yvel(i)     += t_l*vy(j-1, kgrid, rgrid)
        
        dpressure(i) += dt_l*dpres(j-1, kgrid, rgrid)
        dtemp(i)     += dt_l*dtmp(j-1, kgrid,rgrid)
        dpotential(i)+= dt_l*dpot(j-1, kgrid,rgrid)
        dxvel(i)     += dt_l*vx(j-1, kgrid, rgrid)
        dyvel(i)     += dt_l*vy(j-1, kgrid, rgrid)

        d2pressure(i) += d2t_l*dpres(j-1, kgrid, rgrid)
        d2temp(i)     += d2t_l*dtmp(j-1, kgrid, rgrid)
        d2potential(i)+= d2t_l*dpot(j-1, kgrid,rgrid)
        d2xvel(i)     += d2t_l*vx(j-1, kgrid, rgrid)
        d2yvel(i)     += d2t_l*vy(j-1, kgrid, rgrid)

                                ;odd field
        m = n + 1d0
        res = chebyshev_poly(m, zbar)
        t_l  = res[0]
        dt_l = res[1]/zmax
        d2t_l= res[2]/zmax^2
        zvel(i)    +=  t_l*vz(j-1, kgrid, rgrid)
        dzvel(i)   += dt_l*vz(j-1, kgrid, rgrid)
        d2zvel(i)  += d2t_l*vz(j-1, kgrid, rgrid)

     endfor
  endfor 
  
  ;get density from pressure and temperature 
   density  = pressure - temp 
  ddensity  = dpressure - dtemp 
  d2density = d2pressure - d2temp 

  ;error test 
 
;  real*8, parameter :: omegaz=1d0, pi = 2d0*acos(0d0)
;  real*8, parameter :: stefan = 5.67d-5, bigG = 6.67d-8, mstar = 1.99d33, au=1.50d13, kboltz = 1.38d-16, mh=1.66d-24
;  real*8, parameter :: mean_mu = 2.33, rstar = kboltz/mh, sig0=2200d0, T0=120d0, opacity0 = 5d-4, year = 3.17d7
 


 
  s     = dcomplex(growth(kgrid,rgrid), freq(kgrid,rgrid))

  radius = basic2d(0,rgrid)
  Q3d   = basic2d(2,rgrid)

  alpha = alpha_arr(rgrid)
  alpha_b = bulkv*alpha

  dlogrhonu = dlogd ; assume alpha independent of height 
  dlogrhonu_b = dlogrhonu

  print, 'assume opacity scale = 1'
  omega = sqrt(mstar*bigG/(radius*au)^3)
  rhomid = Omega^2d0/(4d0*!dpi*bigG*Q3d)
  Tmid   = temp2d(rgrid)

  fzero  = 16d0*stefan*(gmma-1d0)*Omega
  fzero /=  3d0*opacity0*rhomid^2d0*(rstar/mean_mu)^2d0
  fzero *= Tmid^(2d0 - smallb)

  fzero_thin = 4d0*stefan*(gmma-1d0)*opacity0*temp2d(rgrid)^6d0
  fzero_thin/= (rstar/mean_mu)*temp2d(rgrid)*omega 


  f0    = fzero/(gmma-1d0) 
  bgtest1 = t*dlogt*f0*t^(3d0-smallb)/d/fz + 1d0  
  bgtest2 = t*dlogp + zaxis + sg 
  bgtest3 = dsg - d/Q3d  
  
  ;stop

  
  delta_logrhonu = (1d0+mu)*density + lambda*pressure
  delta_logrhonu/= t
 
  ;density eq 
  test1 = s*density/t + ii*k*xvel + zvel*dlogd + dzvel
  test1/= max(abs(s*density/t))

  ;energy eq 
  
  if(thincool eq -1d0) then begin
;     radiation = -fzero*k*k*(t/d^2)*temp
;     radiation+= fzero*(t/d^2)*( (dlogt - dlogd)*( dtemp + (temp - density)*dlogt ) $
;                                 +d2temp + (dtemp-ddensity)*dlogt + (temp - density)*d2logt)
     radiation = -fzero*k*k*(t^(3d0-smallb)/d^2)*temp
     radiation+=  fzero*(t^(3d0-smallb)/d^2)*( ((3d0-smallb)*dlogt - dlogd)*( dtemp + ((3d0-smallb)*temp - density)*dlogt ) $
                                               + d2temp + ((3d0-smallb)*dtemp-ddensity)*dlogt + ((3d0-smallb)*temp - density)*d2logt)
  endif else begin
     radiation = -(pressure - theta2d(rgrid)*density)/bcool2d(rgrid)
;     radiation -= fzero_thin*(density + 6d0*temp)
  endelse

  dHvisc = alpha*smallq^2*delta_logrhonu - 2d0*ii*k*alpha*smallq*yvel
  dHvisc*= gmma-1d0 

  test2 = s*pressure + zvel*t*dlogp + gmma*t*(ii*k*xvel + dzvel) $
           - dHvisc - radiation 
  test2/= max(abs(s*pressure))

  ;xmom eq 
;  dFx = alpha*( d2xvel + dlogrhonu*dxvel - four_thirds*k^2*xvel )
;  dFx+= ii*k*alpha*(one_third*dzvel + dlogrhonu*zvel)
;  dFx-= alpha_b*k^2*xvel
;  dFx+= ii*k*alpha_b*dzvel 
  
  dFx1 = 2d0*alpha*ii*k*xvel + (alpha_b - 2d0*alpha/3d0)*(ii*k*xvel + dzvel)
  dFx1*= ii*k 
  dFx2 = alpha*dlogrhonu*(ii*k*zvel + dxvel) + alpha*(ii*k*dzvel + d2xvel)
  dFx  = dFx1 + dFx2 
  test3 = s*xvel + ii*k*pressure + ii*k*potential - 2d0*yvel - dFx
  test3/= max(abs(pressure))/zmax

  ;ymom eq 
;  dFy = -ii*k*alpha*smallq*delta_logrhonu
;  dFy+= alpha*( d2yvel + dlogrhonu*dyvel - k^2*yvel)
  dFy1 = -smallq*alpha*delta_logrhonu + alpha*ii*k*yvel 
  dFy1*= ii*k 
  dFy2 = alpha*dlogrhonu*dyvel + alpha*d2yvel
  dFy  = dFy1 + dFy2 
  test4 = s*yvel + (2d0-smallq)*xvel - dFy
  test4/= max(abs(pressure))/zmax

  ;zmom eq 
;  dFz = ii*k*alpha*( one_third*dxvel - two_thirds*dlogrhonu*xvel)
;  dFz+= alpha*( four_thirds*d2zvel + four_thirds*dlogrhonu*dzvel - k^2*zvel)
;  dFz+= alpha_b*(d2zvel + dlogrhonu_b*dzvel)
;  dFz+= ii*alpha_b*k*(dxvel + dlogrhonu_b*xvel)

  dFz1 =ii*k*alpha*(ii*k*zvel + dxvel)
  dFz2 = 2d0*alpha*dlogrhonu*dzvel + 2d0*alpha*d2zvel $
         + (alpha_b - 2d0*alpha/3d0)*(ii*k*dxvel + d2zvel) $
         + (alpha_b*dlogrhonu_b - 2d0*alpha*dlogrhonu/3d0)*(ii*k*xvel + dzvel)
  dFz = dFz1 + dFz2
  test5 = s*zvel + dpressure + dlogd*pressure - density*dlogp + dpotential - dFz
  test5/= max(abs(dpressure)) 

  ;poisson eq
  test6 = d2potential - k*k*potential - d*density/t/Q3d
  test6/= max(abs(d*density/t/Q3d))

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
     bound1/= max(abs(potential))/zmax
  endif
  if(vbc eq 'free') then begin
     bound1 =s*(dpotential(nz-1) + k*potential(nz-1))  + d(nz-1)*zvel(nz-1)/Q3d 
     bound1/= max(abs(s*potential))/zmax
  endif

  bound2 = dxvel(nz-1)
  bound2/= max(abs(xvel))/zmax + 1d-16
  bound3 = dyvel(nz-1)
  bound3/= max(abs(yvel))/zmax + 1d-16

  if(vbc eq 'wall') then begin
     bound4 = zvel(nz-1)
     bound4/= max(abs(zvel)) + 1d-16 
  endif 
  if(vbc eq 'free') then begin
     bound4 = s*pressure(nz-1) + t(nz-1)*dlogp(nz-1)*zvel(nz-1)
     bound4/= max(abs(s*pressure))
  endif 

  ;zero temp bc 
  if(thincool eq -1d0) then begin
     bound5 = temp(nz-1)
     bound5/= max(abs(temp)) + 1d-16 
  endif else bound5 = 0d0 

  boundary_error = max(abs([bound1,bound2,bound3,bound4,bound5]))
  print, 'boundary error', boundary_error

;  density *= dens/csq
  density *= conj(density(0)) 
  density /= max(abs(density)) 
  xtitle = 'z/H'
  ytitle = textoidl('c_s^2\delta\rho/\rho')
  set_plot,'ps'
  file = strcompress('eigenvec_dens.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, zaxis,real_part(density),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=[-1,1],$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval
  oplot, zaxis, imaginary(density), thick=4, linestyle=1 
  device,/close

;  pressure *= dens 
  pressure *= conj(pressure(0))
  pressure /= max(abs(pressure))
  ytitle = textoidl('\deltaP/\rho')
  set_plot,'ps'
  file = strcompress('eigenvec_pres.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, zaxis,real_part(pressure),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=[-1,1],$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval
  oplot, zaxis, imaginary(pressure), thick=4, linestyle=1
  device,/close

  temp *= conj(temp(0))
  temp /= max(abs(temp))
  ytitle = textoidl('\deltaT')
  set_plot,'ps'
  file = strcompress('eigenvec_temp.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, zaxis,real_part(temp),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=[-1,1],$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval
  oplot, zaxis, imaginary(temp), thick=4, linestyle=1
  device,/close

  xvel *= conj(xvel(0))
  xvel /= max(abs(xvel))
  ytitle = textoidl('\deltav_x')
  set_plot,'ps'
  file = strcompress('eigenvec_vx.ps',/remove_all)
  device, filename=file $
          ,bits_per_pixel=8,xsize=8, ysize=4.5,xoffset=0,yoffset=0,/inches,/color
  plot, zaxis,real_part(xvel),xmargin=[8.3,1.7],ymargin=[3.2,1.8], ystyle=0, xstyle=1 $
        ,charsize=1.5, thick=4, xrange=xrange, xtitle=xtitle, yrange=[-1,1],$
        linestyle = 0, ytitle =ytitle, xtickinterval=xtickinterval, ytickinterval=ytickinterval
  oplot, zaxis, imaginary(xvel), thick=4, linestyle=1
  device,/close
end 

