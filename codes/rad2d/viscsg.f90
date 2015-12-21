module global
  implicit none
  character*5 :: method 
  real*8, parameter :: pi = 2d0*acos(0d0)
  real*8, parameter :: stefan = 5.67d-5, bigG = 6.67d-8, mstar = 1.99d33, au=1.50d13, kboltz = 1.38d-16, mh=1.66d-24
  real*8, parameter :: mean_mu = 2.3, rstar = kboltz/mh, opacity_b2 = 5d-4, opacity_b0=0.24, year = 3.16d7 
  real*8 :: fixQ, opscale, k, bcool, gmma, smallq, mu, lambda, alpha, eta, alpha_m, opacity0, opacity_law   
  real*8 :: Htrue, bvisc, vbulk, tirr, theta, fixalpha, bigQ, mdot, radius 
end module global

program viscsg
  use global
  implicit none 
  logical :: gsoft, refine, filter, no_oscil
  integer, parameter :: order = 4
  integer :: nk, nr, i, j, loc(1), m, iflag
  real*8, parameter :: tol=1d-15
  real*8, parameter :: osc_lim = 1d-9
  real*8 :: grav3d, bigQ_fix
  real*8 :: kmin, kmax, rmin, rmax, dk, dlogr, surf, temp  
  real*8 :: c0, c1, c2, c3, c4, growth(order), freq(order), bigA, bigB, bigC, bigD, bigF
  real*8 :: x(2), fvec(2), err(order)
  complex*16 :: coeff(order+1), roots(order), sfreq, det 
  complex*16, allocatable :: rate_fixr(:)
  real*8, allocatable :: raxis(:), kaxis(:), max_rate(:), max_freq(:),  max_k(:), alpha_arr(:), theta_arr(:), Q_arr(:), error(:), error_fixr(:), bcool_arr(:)
  real*8, allocatable :: temp_arr(:), surf_arr(:)  
  real*8, external :: zmag, toomreQ, optical_depth 
  external :: fcn
  namelist /params/ fixQ, mdot, opscale, opacity_law, gsoft, Htrue, gmma, smallq, mu, lambda, fixalpha, alpha_m, vbulk, tirr, refine, filter, no_oscil, method
  namelist /grid/ kmin, kmax, nk, rmin, rmax, nr

  !read input parameters.
  open(7, file="input")
  read(7, nml=params)
  read(7, nml=grid)
  close(7)

  if(opacity_law.eq.0d0) opacity0 = opacity_b0
  if(opacity_law.eq.2d0) opacity0 = opacity_b2


  !convert mdot (sol mass per year) to cgs 
  mdot = mdot*mstar/year 

  if(gsoft .eqv. .true.) then 
     print*, 'soften self-grav due to 3D'
     grav3d = Htrue
  else
     grav3d = 0d0 
  endif

  !output parameters
  open(10,file='params.dat')
  write(10,fmt='(7(e22.15,x))'), fixQ, gmma, mu, lambda, smallq, vbulk, tirr  
  close(10)
  
  !allocate array
  allocate(kaxis(nk))
  allocate(rate_fixr(nk))
  allocate(error_fixr(nk))
  
  allocate(raxis(nr))
  allocate(max_rate(nr))
  allocate(max_freq(nr))
  allocate(max_k(nr))
  allocate(alpha_arr(nr))
  allocate(bcool_arr(nr))
  allocate(Q_arr(nr))
  allocate(theta_arr(nr))
  allocate(error(nr))
  allocate(temp_arr(nr))
  allocate(surf_arr(nr)) 

  if(nk .gt. 1) then 
     dk = (kmax - kmin)/(dble(nk)-1d0)
     do i=1, nk
        kaxis(i) = kmin + (dble(i)-1d0)*dk
     enddo
  else
     kaxis(1) = kmin
  endif
   
  if(nr .gt. 1) then 
     rmin = log10(rmin)
     rmax = log10(rmax)
     dlogr = (rmax-rmin)/(dble(nr)-1d0)
     do i=1, nr
        raxis(i) = 10d0**(rmin + dlogr*(dble(i)-1d0))
     enddo
  else
     raxis(1) = rmin
  endif

  !initial guess for surface density and temp at inner radius 
  surf = 1d3
  temp = 1d3
  
  do i=1, nr     
     radius = raxis(i)
     
     call get_surface_density_and_temperature(surf, temp) 

     surf_arr(i) = surf 
     temp_arr(i) = temp 

     call get_alpha(radius, surf, temp)
     call get_bcool_and_theta(radius, surf, temp) 
     bigQ_fix  = toomreQ(radius, surf, temp)
      
     alpha_arr(i) = alpha  
     theta_arr(i) = theta 
     bcool_arr(i) = bcool 
     Q_arr(i)     = bigQ_fix


    print*, radius, temp, surf !optical_depth(surf,temp)

    if(fixalpha.ge.0d0) alpha = fixalpha  !set alpha in the pert state 

    bvisc        = alpha*vbulk

    eta = 1d0 - alpha*bcool*smallq**2*(gmma-1d0)*lambda
     
     do j=1, nk
        k = kaxis(j)
       
        bigQ = bigQ_fix*( 1d0 + abs(k)*grav3d ) 
  
        c0 = (1d0/bigQ)*alpha*k**2d0*(2d0*eta*(-k + bigQ*(1d0 + mu)*smallq) - &
             alpha*bcool*(-1d0 + gmma)*smallq**2d0*(4d0*k*lambda + &
             bigQ*(1d0 + mu)*(k**2d0 - 2d0*lambda*smallq)) + &
             bigQ*(k**2d0 + 2d0*lambda*smallq)*theta)
        
        c1 = (1d0/(3d0*bigQ))*(-6d0*k*(eta + alpha*bcool*k**2d0) + &
             bigQ*(eta*(12d0 + 4d0*alpha**2d0*k**4d0 + 3d0*alpha*bvisc*k**4d0 - 6d0*smallq) + &
             k**2d0*(8d0*alpha**3d0*bcool*(-1d0 + gmma)*k**2d0*lambda*smallq**2d0 + &
             6d0*alpha**2d0*bcool*bvisc*(-1d0 + gmma)*k**2*lambda*smallq**2d0 + &
             3d0*alpha*bcool*(-(3d0 + mu)*(-2d0 + smallq)*smallq + &
             gmma*(k**2d0 + smallq*(-4d0 + 2d0*lambda + (3d0 + mu)*smallq))) + &
             3d0*theta)))

        c2 = (1d0/(3d0*bigQ))*(bigQ*(7d0*alpha + 3d0*bvisc)*eta*k**2d0 + &
             bcool*(-6d0*k + &
             bigQ*(12d0 + 4d0*alpha**2d0*k**4d0 + 3d0*alpha*bvisc*k**4d0 - 6d0*smallq - &
             6d0*alpha**2d0*k**2d0*lambda*smallq**2d0 + &
             3d0*gmma*k**2d0*(1d0 + 2d0*alpha**2d0*lambda*smallq**2d0))))

        c3 = eta + (1d0/3d0)*bcool*(7d0*alpha + 3d0*bvisc)*k**2d0

        c4 = bcool 

        coeff(1) = dcmplx(c0, 0d0)
        coeff(2) = dcmplx(c1, 0d0)
        coeff(3) = dcmplx(c2, 0d0)
        coeff(4) = dcmplx(c3, 0d0)
        coeff(5) = dcmplx(c4, 0d0)
        
        call get_roots(coeff, roots)

        do m=1, order 
           
           !refine root if desired (using hybrd). iterative solve for poly. roots
           if(refine .eqv. .true.) call refine_roots(roots(m))
          
           !calculate error 
           x(1) = dble(roots(m))
           x(2) = dimag(roots(m))
           
           call fcn(2,x,fvec,iflag)
           det = dcmplx(fvec(1),fvec(2))
           err(m) = zmag(det)
           
           !filter out root with too large error
           if(filter .eqv. .true.) then
              if(err(m) .gt. tol) roots(m) = dcmplx(-1d0,0d0)
           endif
           
           if((no_oscil.eqv..true.).and.(abs(x(2)/x(1)).ge.osc_lim)) then
              roots(m) = dcmplx(-1d0,0d0)
           endif

        enddo
        
        !pick most unstable root  among all roots, and store the associated error
        loc = maxloc(dble(roots)) 
        rate_fixr(j) = roots(loc(1))
        error_fixr(j)= err(loc(1))
     enddo

     !pick most unstable k among all k 
     loc    = maxloc(dble(rate_fixr))
     sfreq =  rate_fixr(loc(1))

     max_rate(i) = dble(sfreq)
     max_freq(i) = dimag(sfreq)
     max_k(i) = kaxis(loc(1))
     error(i) = error_fixr(loc(1))
  enddo
  
  !output
  open(10,file='output.dat')
  do i=1, nr
     write(10,fmt='(11(e22.15,x))'), raxis(i), alpha_arr(i), max_rate(i), max_freq(i), max_k(i), error(i), Q_arr(i), theta_arr(i), bcool_arr(i), surf_arr(i), temp_arr(i)  
  enddo
  close(10)
  
end program viscsg

real*8 function zmag(z)
  implicit none
  complex*16, intent(in) :: z
  real*8 :: re, im
  re = dble(z)
  im = dimag(z)
  zmag = sqrt(re*re + im*im)
end function zmag

real*8 function bigOmega(rad)
  use global
  implicit none
  real*8, intent(in) :: rad

  bigOmega = sqrt(bigG*mstar/au**3d0)*rad**(-3d0/2d0)
end function bigOmega

real*8 function cs2(temp)
use global
implicit none
real*8,intent(in) :: temp 
cs2 = rstar*temp/mean_mu 
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

real*8 function FJ_m(rad, surface_density, temp)
use global 
implicit none 
real*8, intent(in) :: rad, surface_density, temp 
real*8, external :: cs2 

FJ_m = 3d0*pi*alpha_m*cs2(temp)*surface_density*rad*rad*au*au 
end function FJ_m 

real*8 function opacity(temp)
  use global
  implicit none
  real*8, intent(in) :: temp
  
  opacity = opacity0*opscale*temp**opacity_law
end function opacity

real*8 function optical_depth(surface_density, temp)
use global
implicit none
real*8, intent(in) :: temp, surface_density 
real*8, external :: opacity 

optical_depth = opacity(temp)*surface_density 
end function optical_depth 

real*8 function ftau(tau)
implicit none
real*8,intent(in)::tau 
ftau = tau + 1d0/tau 
end function ftau 

real*8 function gtau(tau)
implicit none
real*8,intent(in)::tau
gtau = tau - 1d0/tau
end function gtau


real*8 function toomreQ(rad, surface_density, temp)
  use global
  implicit none
  real*8, intent(in) :: rad, surface_density, temp 
  real*8 :: cs, Omega
  real*8, external :: bigOmega, cs2 

  Omega = bigOmega(rad) 
  cs = sqrt(cs2(temp)) 
  
  toomreQ = cs*Omega/(pi*bigG*surface_density)
end function toomreQ

subroutine get_alpha(rad, surface_density, temp) 
  use global 
  implicit none
  real*8, intent(in) :: rad, surface_density, temp 
  real*8 :: csq, Omega, tau, f_tau, TQ
  real*8, external :: cs2, bigOmega, ftau, optical_depth, T_Q
  
  if(fixalpha .lt. 0d0) then !self-consistent alpha
     csq   = cs2(temp) 
     Omega = bigOmega(rad) 
     tau   = optical_depth(surface_density, temp)
     f_tau = ftau(tau) 
     TQ    = T_Q(radius, surface_density)

     alpha = mdot*Omega/(3d0*pi*csq*surface_density)  
!      alpha = 2d0*stefan*(temp**4d0 - tirr**4d0)/(f_tau*smallq**2d0*Omega*csq*surface_density) 
!    alpha = 8d0*stefan*(pi*bigG*fixQ)**6d0*(mean_mu/rstar)**4d0*surface_density**5d0*(1d0 - tirr**4d0/TQ**4d0)
!    alpha = alpha/(9d0*f_tau*Omega**7d0)
  else !fixed alpha 
     alpha = fixalpha
  endif
end subroutine get_alpha

subroutine get_bcool_and_theta(rad, surface_density, temp)
  use global
  implicit none
  real*8, intent(in) :: rad, surface_density, temp 
  real*8 :: tau, f_tau, g_tau, csq, Omega, c_one, c_two, theta_tilde
  real*8, external :: cs2, optical_depth, ftau, gtau, bigOmega
  
  Omega   = bigOmega(rad) 
  csq     = cs2(temp)

  tau     = optical_depth(surface_density, temp)
  f_tau   = ftau(tau) 
  g_tau   = gtau(tau) 

  theta_tilde = (tirr/temp)**4d0

  c_one = 4d0 - opacity_law*(1d0 - theta_tilde)*g_tau/f_tau
  c_two = 4d0 +(1d0 - opacity_law)*(1d0 - theta_tilde)*g_tau/f_tau

  bcool = f_tau*csq*surface_density*Omega 
  bcool = bcool/( 2d0*(gmma-1d0)*stefan*temp**4d0*c_one )

  theta = c_two/c_one 
end subroutine get_bcool_and_theta

subroutine get_surface_density_and_temperature(surf, temp)
use global
implicit none
real*8, intent(inout) :: surf, temp 
 integer, parameter :: n=2
  integer, parameter :: lwa = n*(3*n+13)
  integer :: info
  real*8, parameter :: tol = 1d-15
  real*8 :: x(n), fvec(n), wa(lwa)
  external :: fcn1

  x(1) = surf
  x(2) = temp 

  call hybrd1(fcn1,n,x,fvec,tol,info,wa,lwa)

  surf = x(1)
  temp = x(2) 
end subroutine get_surface_density_and_temperature

subroutine fcn1(n,x,fvec,iflag)
  use global
  implicit none
  integer :: n,iflag
  real*8 :: x(n),fvec(n)
  real*8 :: surf, temp, Omega, f_tau, tau, fjgt, fjm, norm  
  real*8, external :: bigOmega, optical_depth, ftau, FJ_gt, FJ_m

  surf = x(1)
  temp = x(2) 
  Omega= bigOmega(radius) 
  tau  = optical_depth(surf,temp) 
  f_tau= ftau(tau)   

  fvec(1) = 2d0*stefan*(temp**4d0 - Tirr**4d0)/(smallq**2d0*Omega**2d0*f_tau)/mdot - 1d0/(3d0*pi) 

  fjgt = FJ_gt(radius, surf)
  fjm  = FJ_m(radius, surf, temp)

  norm = (radius*au)**2d0*Omega*mdot 

  fvec(2) = 1d0 - fjm/norm - max(fjgt, 0d0)/norm 

!  print*, fvec(1), fvec(2)


end subroutine fcn1 


subroutine get_roots(coeff, roots)
  use global
  implicit none
  logical, parameter :: polish = .true.
  integer,parameter :: N=4
  integer,parameter :: LWORK = 5*N
  integer :: info 
  complex*16, intent(in) :: coeff(N+1)
  complex*16, intent(out) :: roots(N)
  complex*16 :: eigen(N), VL(N), VR(N), WORK(LWORK)
  complex*16 :: a0, a1, a2, a3, a4, matrix(N,N), one, zero
  real*8 :: RWORK(2*N)

  one = dcmplx(1d0,0d0)
  zero= dcmplx(0d0,0d0)
  
  a0 =  coeff(1)
  a1 = -coeff(2)/a0
  a2 = -coeff(3)/a0
  a3 = -coeff(4)/a0
  a4 = -coeff(5)/a0

  matrix(:,:) = zero
  
  matrix(1,:) = (/a1, a2, a3, a4 /)
  matrix(2,1) = one
  matrix(3,2) = one
  matrix(4,3) = one

  if(method .eq. 'eigen') then
     call zgeev ('N', 'N', N, matrix, N, eigen, VL, N, VR, N, WORK, LWORK, RWORK, INFO)
     roots = 1d0/eigen
  endif
  
  if(method .eq. 'roots') then
     call zroots(coeff,N,roots,polish)
  endif
end subroutine get_roots

subroutine refine_roots(root)
  implicit none
  integer, parameter :: n=2
  integer, parameter :: lwa = n*(3*n+13)
  integer :: info 
  real*8, parameter :: tol = 1d-15
  complex*16, intent(inout) :: root
  real*8 :: re, im, x(n), fvec(n), wa(lwa)
  external :: fcn 

  re = dble(root)
  im = dimag(root)
  
  x(1) = re
  x(2) = im 
  
  call hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)

end subroutine refine_roots

subroutine fcn(n,x,fvec,iflag)
  use global
  implicit none 
  integer :: n,iflag
  real*8 :: x(n),fvec(n)
  complex*16 :: sfreq, bigF, bigA, bigB, bigC, bigD, det, bigE
  
  sfreq = dcmplx(x(1),x(2))
  
  bigF = k**2d0*(gmma-1d0)*bcool
  bigF = bigF/(eta + bcool*sfreq)
           
  bigE = gmma/(gmma-1d0) + alpha*smallq**2d0*(1d0+mu)/sfreq + theta/bcool/sfreq/(gmma-1d0)


  bigA =((4d0/3d0)*alpha + bvisc)*k**2d0 + sfreq + bigF*bigE - 2d0*abs(k)/bigQ/sfreq
           
  bigB = 2d0*(alpha*smallq*bigF - 1d0)
           
  bigC = (2d0-smallq) + alpha*smallq*k**2d0*(1d0+mu)/sfreq + alpha*smallq*lambda*bigF*bigE
           
  bigD = alpha*k**2d0 + sfreq + 2d0*alpha**2d0*smallq**2d0*lambda*bigF 
           
  det = bigA*bigD - bigB*bigC 
        
  fvec(1) = dble(det)
  fvec(2) = dimag(det)
  
  return
end subroutine fcn

SUBROUTINE laguer(a,m,x,its)
  implicit none
  INTEGER :: m,its,MAXIT,MR,MT
  REAL*8 :: EPSS
  COMPLEX*16 :: a(m+1),x
  PARAMETER (EPSS=1d-15,MR=8,MT=500,MAXIT=MT*MR)
  INTEGER :: iter,j
  REAL*8 :: abx,abp,abm,err,frac(MR)
  COMPLEX*16 :: dx,x1,b,d,f,g,h,sq,gp,gm,g2
  SAVE frac
  DATA frac /.5d0,.25d0,.75d0,.13d0,.38d0,.62d0,.88d0,1d0/
  real*8, external :: zmag 
  do iter=1,MAXIT
     its=iter
     b=a(m+1)
     err=zmag(b)
     d=dcmplx(0d0,0d0)
     f=dcmplx(0d0,0d0)
     abx=zmag(x)
     do j=m,1,-1
        f=x*f+d
        d=x*d+b
        b=x*b+a(j)
        err=zmag(b)+abx*err
     enddo
     err=EPSS*err
     if(zmag(b).le.err) then
        return
     else
        g=d/b
        g2=g*g
        h=g2-2d0*f/b
        sq=sqrt((m-1)*(m*h-g2))
        gp=g+sq
        gm=g-sq
        abp=abs(gp)
        abm=abs(gm) 
        if(abp.lt.abm) gp=gm
        if (max(abp,abm).gt.0d0) then
           dx=m/gp
        else
           dx=exp(dcmplx(log(1d0+abx),dble(iter)))
        endif
     endif
     x1=x-dx
     if(x.eq.x1)return
     if (mod(iter,MT).ne.0) then
        x=x1
     else
        x=x-dx*frac(iter/MT)
     endif
  enddo
  pause 'too many iterations in laguer'
  return
END SUBROUTINE laguer

SUBROUTINE zroots(a,m,roots,polish)
  implicit none 
  INTEGER :: m,MAXM
  REAL*8 :: EPS
  COMPLEX*16 :: a(m+1),roots(m)
  LOGICAL :: polish
  PARAMETER (EPS=1d-15,MAXM=101)
  INTEGER :: i,j,jj,its
  COMPLEX*16 :: ad(MAXM),x,b,c
  do j=1,m+1
     ad(j)=a(j)
  enddo
  do j=m,1,-1
     x=dcmplx(0d0,0d0)
     call laguer(ad,j,x,its)
     if(abs(aimag(x)).le.2d0*EPS**2*abs(dble(x))) x=dcmplx(dble(x),0d0)
     roots(j)=x
     b=ad(j+1)
     do  jj=j,1,-1
        c=ad(jj)
        ad(jj)=b
        b=x*b+c
     enddo
  enddo
  if (polish) then
     do j=1,m
        call laguer(a,m,roots(j),its)
     enddo
  endif
  do j=2,m
     x=roots(j)
     do i=j-1,1,-1
        if(dble(roots(i)).le.dble(x))goto 10
        roots(i+1)=roots(i)
     enddo
     i=0
10   roots(i+1)=x
  enddo
  return
END SUBROUTINE zroots

