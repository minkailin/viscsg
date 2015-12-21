subroutine get_thickness_and_Q3d
  use global
  implicit none
  integer, parameter :: n = 2
  integer, parameter :: lwa = n*(3*n+13)
  integer            :: info 
  real*8,  parameter :: guess = 1d0, tol=1d-15
  real*8             :: x(n), fvec(n), wa(lwa)
  external           :: fcn 
  
  x(1) = guess
  x(2) = Q2d/4d0

  call hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
  
  bigH = x(1)
  Q3d  = x(2) 
end subroutine get_thickness_and_Q3d

subroutine fcn(n,x,fvec,iflag)
  use global
  implicit none
  integer, parameter :: neq = 4, itol = 1, itask=1, iopt=0, MF = 10 
  integer, parameter :: lrw =2*(20 + 16*NEQ), liw = 60
  integer            :: n, iflag, istate, iwork(liw), IPAR
  real*8, parameter  :: rtol = 1d-15, atol = 1d-15
  real*8             :: x(n), fvec(n)
  real*8             :: Y(neq), Z, ZOUT, rwork(lrw), RPAR 
  real*8             :: final_density
  external           :: F, JAC
 
  Z    = 0d0 !start from midplane
  ZOUT = x(1)!guess end point 
  Q3d  = x(2)!3D Q parameter  

  Y(1) = 1d0 !normalized pressure at z=0
  Y(2) = 1d0 !normalized temperature at z=0
  Y(3) = 0d0 !flux at the midplane
  Y(4) = 0d0 !self-gravity force at midplane 
  
  istate = 1 
  call DVODE (F, NEQ, Y, Z, ZOUT, ITOL, RTOL, ATOL, ITASK, &
       ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF, &
       RPAR, IPAR)

  !test final density and Q2d 
  final_density = Y(1)/Y(2)
  fvec(1) = final_density - (1d0 - zmax)
  fvec(2) = Y(4) - 2d0/Q2d 
end subroutine fcn

subroutine F (NEQ, Z, Y, YDOT, RPAR, IPAR)
  use global
  implicit none
  integer :: NEQ, IPAR
  real*8  :: Z, Y(NEQ), YDOT(NEQ), RPAR
  real*8  :: pressure, rho, temperature, g, fz
  real*8  :: rho_mid, Omega, F0
  real*8, external :: bigOmega 

  pressure     = Y(1)
  temperature  = Y(2)
  fz           = Y(3)
  g            = Y(4) 
  rho          = pressure/temperature

  Omega   = bigOmega(radius) 
  rho_mid = Omega**2d0/(4d0*pi*bigG*Q3d)
  
  F0  = 3d0*opacity0*opscale*rho_mid**2d0*(rstar/mean_mu)**2d0
  F0  = F0/(16d0*stefan*Omega*Tmid**(2d0-opacity_law))

  YDOT(1) = -rho*(omegaz**2d0*Z + g)
  YDOT(2) = -F0*rho*fz/temperature**(3d0-opacity_law)
  YDOT(3) = alpha*smallq**2d0*rho
  YDOT(4) = rho/Q3d 

end subroutine F

SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
  DOUBLE PRECISION T, Y(NEQ), PD(NROWPD,NEQ), RPAR
end SUBROUTINE JAC

subroutine get_vertical_structure(zout, d, dlogd, p, dlogp, tmp, dlogtmp, d2logtmp, grav, dgrav, fz)
  use global
  implicit none
  real*8, intent(in) :: zout
  real*8, intent(out):: d, dlogd, p, dlogp, tmp, dlogtmp, d2logtmp, grav, dgrav, fz 
  integer, parameter :: neq = 4, itol = 1, itask=1, iopt=0, MF = 10 
  integer, parameter :: lrw =2*(20 + 16*NEQ), liw = 60
  integer            :: istate, iwork(liw), IPAR
  real*8, parameter  :: rtol = 1d-15, atol = 1d-15
  real*8             :: Y(neq), Z, rwork(lrw), RPAR, YDOT(neq), F0, Omega, rho_mid 
  real*8, external   :: bigOmega 
  external           :: F, JAC
  
  Omega   = bigOmega(radius) 
  rho_mid = Omega**2d0/(4d0*pi*bigG*Q3d)
  
  F0  = 3d0*opacity0*opscale*rho_mid**2d0*(rstar/mean_mu)**2d0
  F0  = F0/(16d0*stefan*Omega*Tmid**(2d0-opacity_law))

  Z    = 0d0 !start from midplane
  
  Y(1) = 1d0 !normalized pressure at z=0
  Y(2) = 1d0 !normalized temperature at z=0
  Y(3) = 0d0 !flux at the midplane
  Y(4) = 0d0 !self-gravity force at midplane 

  istate = 1 
  call DVODE (F, NEQ, Y, Z, ZOUT, ITOL, RTOL, ATOL, ITASK, &
       ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF, &
       RPAR, IPAR)

  !assign density and get log derivative
  p   = Y(1)
  tmp = Y(2) 
  d   = p/tmp 
  fz  = Y(3)
  grav= Y(4) 

  call F (NEQ, ZOUT, Y, YDOT, RPAR, IPAR)
  dlogp    = YDOT(1)/p
  dlogtmp  = YDOT(2)/tmp
  dlogd    = dlogp - dlogtmp  
  d2logtmp = (dlogd*fz - (4d0-opacity_law)*dlogtmp*fz + YDOT(3))*(-F0*d/tmp**(4d0-opacity_law))
  dgrav    = YDOT(4) 
end subroutine get_vertical_structure

subroutine basic_setup(ncells)
  use global
  implicit none
  character*3 :: nbstring
  integer, intent(in) :: ncells
  integer :: i , j, jmid
  real*8 :: dz, zout, d, dlogd, p, dlogp, tmp, dlogtmp, d2logtmp, grav, dgrav, fz
  real*8, external :: gviscosity

  !allocate arrays for basic state 
  allocate(zaxis(ncells))
  allocate(dens(ncells))
  allocate(dlogdens(ncells))
  allocate(pres(ncells))
  allocate(dlogpres(ncells))
  allocate(temp(ncells))
  allocate(dlogtemp(ncells))
  allocate(d2logtemp(ncells))
  allocate(sg(ncells))
  allocate(dsg(ncells))
  allocate(radflux(ncells))

  if(ncells.eq.nz) then !set up basic state on spectral grid
     lmax = 2*Nz - 1 !maximum chebyshev order (odd)
     jmid = Nz + 1
     do i = 1, nz
        j = i + jmid - 1
        zaxis(i) = -bigH*cos(pi*(j-1d0)/lmax)
     enddo
  endif
  
  if(ncells.eq.nzout) then !setup basic state on fine grid for output 
     dz = bigH/(dble(nzout)-1d0)
     do i = 1, nzout
        zaxis(i) = (dble(i)-1d0)*dz
     enddo
  endif
  
  !setup vertical structure grids 
  do i=1, ncells 
     zout = zaxis(i)
     call get_vertical_structure(zout, d, dlogd, p, dlogp, tmp, dlogtmp, d2logtmp, grav, dgrav, fz)
     dens(i)     = d
     dlogdens(i) = dlogd
     pres(i)     = p
     dlogpres(i) = dlogp 
     temp(i)     = tmp
     dlogtemp(i) = dlogtmp  
     d2logtemp(i)= d2logtmp 
     sg(i)       = grav
     dsg(i)      = dgrav 
     radflux(i)  = fz
  enddo
   
end subroutine basic_setup
