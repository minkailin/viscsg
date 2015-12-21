real*8 function bigOmega(rad)
  use global
  implicit none
  real*8, intent(in) :: rad

  bigOmega = sqrt(bigG*mstar/au**3d0)*rad**(-3d0/2d0)
end function bigOmega

real*8 function cs2(tmp)
  use global
  implicit none
  real*8,intent(in) :: tmp

  cs2 = rstar*tmp/mean_mu
end function cs2

real*8 function T_Q(rad, surface_density)
  use global
  implicit none
  real*8, intent(in) :: rad, surface_density
  real*8 :: Omega
  real*8, external :: bigOmega

  Omega = bigOmega(rad)
  T_Q = (mean_mu/rstar)*(pi*bigG*surface_density*fixQ/Omega)**2d0
end function T_Q

real*8 function FJ_gt(rad, surface_density)
  use global
  implicit none
  real*8, intent(in) :: rad, surface_density
  real*8 :: Omega, tauQ, ftauQ, TQ
  real*8, external :: bigOmega, T_Q, optical_depth, ftau

  Omega = bigOmega(rad)
  TQ   = T_Q(rad, surface_density)
  tauQ = optical_depth(surface_density,TQ)
  ftauQ= ftau(tauQ)

  FJ_gt = 8d0*pi*rad*rad*au*au/(3d0*Omega)
  FJ_gt = FJ_gt*stefan*(TQ**4d0-Tirr**4d0)/ftauQ
end function FJ_gt

real*8 function FJ_m(rad, surface_density, tmp)
  use global
  implicit none
  real*8, intent(in) :: rad, surface_density, tmp
  real*8, external :: cs2
  
  FJ_m = 3d0*pi*alpha_m*cs2(tmp)*surface_density*rad*rad*au*au
end function FJ_m

real*8 function opacity(tmp)
  use global
  implicit none
  real*8, intent(in) :: tmp
  
  opacity = opacity0*opscale*tmp**opacity_law
end function opacity

real*8 function optical_depth(surface_density, tmp)
  use global
  implicit none
  real*8, intent(in) :: tmp, surface_density
  real*8, external :: opacity
  optical_depth = opacity(tmp)*surface_density
end function optical_depth

real*8 function ftau(tau)
  implicit none
  real*8,intent(in)::tau
  ftau = tau + 1d0/tau
end function ftau

real*8 function toomreQ(rad, surface_density, tmp)
  use global
  implicit none
  real*8, intent(in) :: rad, surface_density, tmp
  real*8 :: cs, Omega
  real*8, external :: bigOmega, cs2
  
  Omega = bigOmega(rad)
  cs = sqrt(cs2(tmp))
  
  toomreQ = cs*Omega/(pi*bigG*surface_density)
end function toomreQ

subroutine get_alpha(rad, surface_density, tmp)
  use global
  implicit none
  real*8, intent(in) :: rad, surface_density, tmp
  real*8 :: csq, Omega, tau, f_tau
  real*8, external :: cs2, bigOmega, ftau, optical_depth
  
  csq   = cs2(tmp)
  Omega = bigOmega(rad)
  tau   = optical_depth(surface_density, tmp)
  f_tau = ftau(tau)
  
  alpha = mdot*Omega/(3d0*pi*csq*surface_density)  
  !alpha = 2d0*stefan*(tmp**4d0 - tirr**4d0)/(f_tau*smallq**2d0*Omega*csq*surface_density)
end subroutine get_alpha

subroutine get_surface_density_and_temperature(surf, tmp)
  use global
  implicit none
  real*8, intent(inout) :: surf, tmp
  integer, parameter :: n=2
  integer, parameter :: lwa = n*(3*n+13)
  integer :: info
  real*8, parameter :: tol = 1d-9
  real*8 :: x(n), fvec(n), wa(lwa)
  external :: fcn1
  
  x(1) = surf
  x(2) = tmp
  
  call hybrd1(fcn1,n,x,fvec,tol,info,wa,lwa)
  
  surf  = x(1)
  tmp   = x(2)
end subroutine get_surface_density_and_temperature

subroutine fcn1(n,x,fvec,iflag)
  use global
  implicit none
  integer :: n,iflag
  real*8 :: x(n),fvec(n)
  real*8 :: surf, tmp, Omega, f_tau, tau, fjgt, fjm, norm
  real*8, external :: bigOmega, optical_depth, ftau, FJ_gt, FJ_m
  
  surf  = x(1)
  tmp   = x(2)
  Omega= bigOmega(radius)
  tau  = optical_depth(surf,tmp)
  f_tau= ftau(tau)
  
  fvec(1) = 2d0*stefan*(tmp**4d0 - Tirr**4d0)/(smallq**2d0*Omega**2d0*f_tau)/mdot - 1d0/(3d0*pi)
  
  fjgt = FJ_gt(radius, surf)
  fjm  = FJ_m(radius, surf, tmp)
  
  norm = (radius*au)**2d0*Omega*mdot
  
  fvec(2) = 1d0 - fjm/norm - max(fjgt, 0d0)/norm
  
end subroutine fcn1

real*8 function gtau(tau)
  implicit none
  real*8,intent(in)::tau
  gtau = tau - 1d0/tau
end function gtau

subroutine get_bcool_and_theta(rad, surface_density, tmp, bcool, theta)
  use global
  implicit none
  real*8, intent(in) :: rad, surface_density, tmp
  real*8, intent(out) :: bcool, theta 
  real*8 :: tau, f_tau, g_tau, csq, Omega, c_one, c_two, theta_tilde
  real*8, external :: cs2, optical_depth, ftau, gtau, bigOmega

  Omega   = bigOmega(rad)
  csq     = cs2(tmp)

  tau     = optical_depth(surface_density, tmp)
  f_tau   = ftau(tau)
  g_tau   = gtau(tau)

  theta_tilde = (tirr/tmp)**4d0

  c_one = 4d0 - 2d0*(1d0 - theta_tilde)*g_tau/f_tau
  c_two = 4d0 -     (1d0 - theta_tilde)*g_tau/f_tau

  bcool = f_tau*csq*surface_density*Omega
  bcool = bcool/( 2d0*(gmma-1d0)*stefan*tmp**4d0*c_one )

  theta = c_two/c_one
end subroutine get_bcool_and_theta
