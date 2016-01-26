module global
  implicit none
  character*5 :: method 
  real*8 :: bigQ, k, bcool, gmma, smallq, mu, lambda, alpha, eta 
  real*8 :: Htrue, bvisc, vbulk, tirr, theta, bigQ_marg2d, bigQ_marg3d
end module global

program viscsg
  use global
  implicit none 
  logical :: gsoft, refine, filter, gvisc 
  integer, parameter :: order = 4
  integer :: nk, nb, i, j, loc(1), m, iflag
  real*8 :: tol = 1d-15
  real*8 :: bigQ_fix, fixalpha, alim, grav3d 
  real*8 :: kmin, kmax, bmin, bmax, dk, dlogb, btrans
  real*8 :: c0, c1, c2, c3, c4, growth(order), freq(order), bigA, bigB, bigC, bigD, bigF
  real*8 :: x(2), fvec(2), err(order)
  complex*16 :: coeff(order+1), roots(order), sfreq, det 
  complex*16, allocatable :: rate_fixb(:)
  real*8, allocatable :: baxis(:), kaxis(:), max_rate(:), max_freq(:),  max_k(:), alpha_arr(:), error(:), error_fixb(:), Q_arr(:) 
  real*8, external :: zmag, gviscosity 
  external :: fcn 
  namelist /params/ bigQ_fix, gsoft, Htrue, gmma, smallq, mu, lambda, fixalpha, gvisc, vbulk, tirr, refine, filter, method
  namelist /grid/ kmin, kmax, nk, bmin, bmax, nb

  !read input parameters.
  open(7, file="input")
  read(7, nml=params)
  read(7, nml=grid)
  close(7)

  if(gsoft .eqv. .true.) then 
     print*, 'soften self-grav due to 3D'
     grav3d = Htrue
  else
     grav3d = 0d0 
  endif

  !input bigQ_fix is initial Q in units of the marginally stable value
     k = -(1d0/6d0) + gmma/( &
          3d0*2d0**(2d0/3d0)*(432d0*gmma**2 - 2d0*gmma**3 - 216d0*gmma**2*smallq + &
          sqrt(-4d0*gmma**6 + (432d0*gmma**2 - 2d0*gmma**3 - 216d0*gmma**2*smallq)**2))**( &
          1d0/3d0)) + (432d0*gmma**2 - 2d0*gmma**3 - 216d0*gmma**2d0*smallq + &
          sqrt(-4d0*gmma**6 + (432d0*gmma**2 - 2d0*gmma**3 - 216d0*gmma**2*smallq)**2))**(1d0/3d0)/( &
          6d0*2d0**(1d0/3d0)*gmma)
     bigQ_marg3d = 1d0/gmma/k/(1d0+k)**2d0

     bigQ_marg2d = 1d0/sqrt(2d0*gmma*(2d0-smallq))

  if(gsoft .eqv. .true.) then
     bigQ_fix = bigQ_fix*bigQ_marg3d
  else 
     bigQ_fix = bigQ_fix*bigQ_marg2d
  endif

  !output parameters
  open(10,file='params.dat')
  write(10,fmt='(7(e22.15,x))'), bigQ_fix, gmma, mu, lambda, smallq, vbulk, tirr  
  close(10)
  
  !allocate array
  allocate(kaxis(nk))
  allocate(baxis(nb))
  allocate(rate_fixb(nk))
  allocate(error_fixb(nk))
  allocate(max_rate(nb))
  allocate(max_freq(nb))
  allocate(max_k(nb))
  allocate(alpha_arr(nb))
  allocate(error(nb))
  allocate(Q_arr(nb))

  if(nk .gt. 1) then 
     dk = (kmax - kmin)/(dble(nk)-1d0)
     do i=1, nk
        kaxis(i) = kmin + (dble(i)-1d0)*dk
     enddo
  else
     kaxis(1) = kmin
  endif
  
!  if(gvisc.eqv..true.) then 
!     if(fixalpha.lt.0d0) then !re-set bmin because we need alpha < 1 
!        bmin = 1d0 - tirr
!        bmin = bmin/( (gmma-1d0)*smallq**2d0 )
!     endif
!  endif
  
  if(nb .gt. 1) then 
     bmin = log10(bmin)
     bmax = log10(bmax)
     dlogb = (bmax-bmin)/(dble(nb)-1d0)
     do i=1, nb
        baxis(i) = 10d0**(bmin + dlogb*(dble(i)-1d0))
     enddo
  else
     baxis(1) = bmin
  endif
  
  do i=1, nb
     bcool = baxis(i)
     
     if(fixalpha .ge. 0d0) then
        alpha = fixalpha
        theta = tirr 
     else 
        theta = tirr
        alpha = 1d0 - theta
        alpha =alpha/( (gmma-1d0)*bcool*smallq**2 )
     endif
     bvisc        = alpha*vbulk
     alpha_arr(i) = alpha 

     if(gvisc.eqv..true.) then !set Q according to alpha 
        bigQ_fix = gviscosity(alpha)
        !print*, 'reset bigQ_fix to', bigQ_fix
     endif
     Q_arr(i)     = bigQ_fix

     eta = 1d0 - alpha*bcool*smallq**2*(gmma-1d0)*lambda
     
     rate_fixb(:) = 0d0  
     
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
!           if(filter .eqv. .true.) then
!              if(err(m) .gt. tol) roots(m) = dcmplx(-1d0,0d0)
!           endif
             
            !filter out over-stable modes 
           if(filter .eqv. .true.) then
              if(abs(x(2)/x(1)) .gt. tol) roots(m) = dcmplx(-1d0,0d0)
           endif


        enddo
        
        !pick most unstable root  among all roots, and store the associated error
        loc = maxloc(dble(roots)) 
        rate_fixb(j) = roots(loc(1))
        error_fixb(j)= err(loc(1))
     enddo

     !pick most unstable k among all k 
     loc    = maxloc(dble(rate_fixb))
     sfreq =  rate_fixb(loc(1))

     max_rate(i) = dble(sfreq)
     max_freq(i) = dimag(sfreq)
     max_k(i) = kaxis(loc(1))
     error(i) = error_fixb(loc(1))
  enddo
  
  !output
  open(10,file='output.dat')
  do i=1, nb
     write(10,fmt='(7(e22.15,x))'), baxis(i), alpha_arr(i), max_rate(i), max_freq(i), max_k(i), error(i), Q_arr(i)
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

real*8 function gviscosity(avisc) !return Q given viscosity
  use global
  implicit none
  real*8, parameter :: amin = 1d-4 
  real*8, intent(in) :: avisc
  
  gviscosity = bigQ_marg2d/sqrt(avisc)
end function gviscosity
  
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
  
!  print*, 'det', det


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

