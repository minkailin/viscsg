module global
  logical :: maximize, no_decay, no_oscil, gvisc 
  character*4 :: vbc 
  integer, parameter :: nvar = 5
  integer :: nz, bignz, nb, nk, lmax, nzout 
  real*8, parameter :: omegaz=1d0, pi = 2d0*acos(0d0)
  real*8 :: Q3d, Q2d, smallq, gmma, bgmma, zmax, bigH, mu, lambda, bmin, bmax, kmin, kmax  
  real*8 :: fixalpha, tirr, bulkv 
  real*8, allocatable :: bvals(:), kvals(:) 
  real*8, allocatable :: zaxis(:), dens(:), dlogdens(:), csq(:), dcsq(:)
  real*8, allocatable :: alpha2d(:,:), dalpha2d(:,:), Qarr(:) 
  real*8, allocatable :: T(:,:), Tp(:,:), Tpp(:,:)
  real*8, allocatable :: T_odd(:,:), Tp_odd(:,:), Tpp_odd(:,:)
  real*8, allocatable :: growth(:,:), growth_fixb(:), freq(:,:)
  complex*16, parameter :: ii = (0d0, 1d0) 
  complex*16, allocatable :: drho(:), dpres(:), dpot(:), vx(:), vy(:), vz(:), invT(:,:)
  complex*16, allocatable :: bigLmatrix(:,:), bigRmatrix(:,:), eigen_vec(:)
  complex*16, allocatable :: L11(:,:), L12(:,:), L13(:,:), L14(:,:), L15(:,:) 
  complex*16, allocatable :: L21(:,:), L22(:,:), L23(:,:), L24(:,:), L25(:,:) 
  complex*16, allocatable :: L31(:,:), L32(:,:), L33(:,:), L34(:,:), L35(:,:) 
  complex*16, allocatable :: L41(:,:), L42(:,:), L43(:,:), L44(:,:), L45(:,:) 
  complex*16, allocatable :: L51(:,:), L52(:,:), L53(:,:), L54(:,:), L55(:,:) 
  complex*16, allocatable :: R11(:,:), R12(:,:), R13(:,:), R14(:,:), R15(:,:) 
  complex*16, allocatable :: R21(:,:), R22(:,:), R23(:,:), R24(:,:), R25(:,:)
  complex*16, allocatable :: R31(:,:), R32(:,:), R33(:,:), R34(:,:), R35(:,:)
  complex*16, allocatable :: R41(:,:), R42(:,:), R43(:,:), R44(:,:), R45(:,:) 
  complex*16, allocatable :: R51(:,:), R52(:,:), R53(:,:), R54(:,:), R55(:,:) 
end module global

program viscsg
  use global
  implicit none
  integer :: i,j, jmid, kcount, bcount, info, lwork
  integer :: loc(1) 
  integer, allocatable :: ipiv(:)
  real*8 :: dk, dlogb, optimize 
  real*8 :: k, zbar, m, T_l, dT_l, d2T_l, n 
  complex*16, allocatable :: work(:)
  namelist /params/ Q3d, smallq, bgmma, gmma, bmin, bmax, nb, kmin, kmax, nk, mu, lambda, gvisc, fixalpha, bulkv, tirr
  namelist /grid/   nz, zmax, nzout, vbc 
  namelist /job/    maximize, no_decay, no_oscil 
  
  !read parameters
  open(7, file="input")
  read(7, nml=params)
  read(7, nml=grid)
  read(7, nml=job)
  close(7)

  !are we maximizing growth rates over k ?
  if(maximize.eqv..true.) then
    optimize = 1d0
  else
    optimize= -1d0
  endif      

  !output parameters for later analysis
  open(10,file='params.dat')
  write(10,fmt='(11(e22.15,x))') Q3d, smallq, bgmma, gmma, mu, lambda, tirr, fixalpha, dble(nz), optimize, bulkv 
  close(10)
  
  !make sure grid cells are odd 
  if(mod(nz,2).eq.0) then
     print*, 'Nz needs to be odd but Nz=', nz
     stop
  endif
  
  !set up parameter survey axis 
  allocate(bvals(nb))
  allocate(kvals(nk))

  !wavenumber space
  if(nk .gt. 1) then
!     dk = (kmax - kmin)/(dble(nk)-1d0)
      dk = log10(kmax/kmin)/(dble(nk)-1d0)
     do i=1, nk 
!        kvals(i) = kmin + dk*(dble(i)-1d0)
         kvals(i) = 10d0**(log10(kmin)+dk*(dble(i)-1d0))
     enddo
  else
     kvals(1) = kmin
  endif
  open(10,file='kvals.dat')
  do i=1, nk
     write(10,fmt='(e22.15)') kvals(i)
  enddo
  close(10)
  
  !cooling space
  if(nb .gt. 1) then
     dlogb = log10(bmax/bmin)/(dble(nb)-1d0)
     do i=1, nb
        bvals(i) = 10d0**(log10(bmin) + dlogb*(dble(i)-1d0))
     enddo
  else
     bvals(1) = bmin 
  endif

  open(10,file='bvals.dat')
  do i=1, nb
     write(10,fmt='(e22.15)') bvals(i)
  enddo
  close(10)
  
  !****************************************!
  !***********BASIC STATE SETUP************!
  !****************************************!
  write(6,fmt='(A)') '************BASIC STATE SETUP***********'
  
  !figure out the vertical domain size and set up vertical grid (z > 0)
  call get_thickness
  write(6,fmt='(A,e22.15)') 'bigH=', bigH 

 !setup basic state and output (fine grid for data analysis)
  call basic_setup(nzout)

  !setup basic state again (coarse spectral grid for eigen problem) 
  call basic_setup(nz) 
  Q2d = -2d0/(csq(nz)*dlogdens(nz) + omegaz**2d0*bigH)
  write(6,fmt='(A,e22.15)') 'Q_2d=', Q2d

  !****************************************!
  !***********EIGEN PROBLEM SETUP**********!
  !****************************************!
  !spectral basis, even for W, phi, vx, vy 
  !                odd for vz                 
  allocate(T(nz,nz))
  allocate(Tp(nz,nz))
  allocate(Tpp(nz,nz))

  allocate(T_odd(nz,nz))
  allocate(Tp_odd(nz,nz))
  allocate(Tpp_odd(nz,nz))
  
  do j=1, nz !jth physical grid
     zbar = zaxis(j)/bigH
     do i=1, nz !ith basis

        !fill the even basis 
        n = 2d0*dble(i-1)
        call chebyshev_poly(n, zbar, T_l, dT_l, d2T_l)
        T(j,i)  =   T_l
        Tp(j,i) =  dT_l/bigH             !convert to deriv wrt physical grid
        Tpp(j,i)= d2T_l/bigH**2d0        !convert to deriv wrt physical grid  
                
        !fill the odd basis (these operate on vz)
        m = n+1d0
        call chebyshev_poly(m, zbar, T_l, dT_l, d2T_l)
        T_odd(j,i)  =   T_l
        Tp_odd(j,i) =  dT_l/bigH
        Tpp_odd(j,i)= d2T_l/bigH**2d0
        
     enddo
  enddo
  
  !we need the inverse of T later (declare as complex)  
  lwork = nz 
  allocate(invT(nz,nz))
  allocate(ipiv(nz))
  allocate(work(lwork))
  invT = T 
  call zgetrf(nz, nz, invT, nz, IPIV, INFO)
  call ZGETRI(nz, invT, nz, IPIV, WORK, LWORK, INFO)

  !for the energy/pressure equation
  allocate(L11(nz,nz)) !operators on W
  allocate(L12(nz,nz)) !operators on phi
  allocate(L13(nz,nz)) !operators on vx
  allocate(L14(nz,nz)) !operators on vy
  allocate(L15(nz,nz)) !operators on vz
  
  allocate(R11(nz,nz)) ; allocate(R12(nz,nz)) ; allocate(R13(nz,nz)) ; allocate(R14(nz,nz)) ; allocate(R15(nz,nz))

  !for density/potential equation
  allocate(L21(nz,nz)) ; allocate(L22(nz,nz)) ; allocate(L23(nz,nz)) ; allocate(L24(nz,nz)) ; allocate(L25(nz,nz))
  allocate(R21(nz,nz)) ; allocate(R22(nz,nz)) ; allocate(R23(nz,nz)) ; allocate(R24(nz,nz)) ; allocate(R25(nz,nz))

  !for the x momentum equation
  allocate(L31(nz,nz)) ; allocate(L32(nz,nz)) ; allocate(L33(nz,nz)) ; allocate(L34(nz,nz)) ; allocate(L35(nz,nz))
  allocate(R31(nz,nz)) ; allocate(R32(nz,nz)) ; allocate(R33(nz,nz)) ; allocate(R34(nz,nz)) ; allocate(R35(nz,nz))

  !for the y momentum equation
  allocate(L41(nz,nz)) ; allocate(L42(nz,nz)) ; allocate(L43(nz,nz)) ; allocate(L44(nz,nz)) ; allocate(L45(nz,nz))
  allocate(R41(nz,nz)) ; allocate(R42(nz,nz)) ; allocate(R43(nz,nz)) ; allocate(R44(nz,nz)) ; allocate(R45(nz,nz))

  !for the z momentum equation
  allocate(L51(nz,nz)) ; allocate(L52(nz,nz)) ; allocate(L53(nz,nz)) ; allocate(L54(nz,nz)) ; allocate(L55(nz,nz))
  allocate(R51(nz,nz)) ; allocate(R52(nz,nz)) ; allocate(R53(nz,nz)) ; allocate(R54(nz,nz)) ; allocate(R55(nz,nz))

  bignz = nvar*nz
  allocate(bigLmatrix(bignz,bignz))
  allocate(bigRmatrix(bignz,bignz))
  allocate(growth(nk,nb))
  allocate(growth_fixb(nk))
  allocate(freq(nk,nb))
  allocate(eigen_vec(bignz))
  allocate(drho(nz))
  allocate(dpres(nz))
  allocate(dpot(nz))
  allocate(vx(nz))
  allocate(vy(nz))
  allocate(vz(nz))

  write(6,fmt='(A)') '************EIGENVALUE PROBLEM**********'
  !loop over parameter space and do eigenvalue problem 
  open(10,file='eigenvalues.dat')
  open(20,file='eigenvectors.dat')
  do bcount=1, nb
     write(6,fmt='(A,e8.1,A)') '************BCOOL=', bvals(bcount),'**************'
     do kcount=1, nk
        call fill_matrices(kcount,bcount)
        call eigenvalue_problem(kcount,bcount)
        
        growth_fixb(kcount) = growth(kcount,bcount) 

        if(maximize.eqv..false.) then 
        write(10,fmt='(3(e22.15,x))') growth(kcount,bcount), freq(kcount,bcount), kvals(kcount) 
        
        dpres(:) = eigen_vec(1:nz)
        dpot(:)  = eigen_vec(nz+1:2*nz)
        vx(:)    = eigen_vec(2*nz+1:3*nz)
        vy(:)    = eigen_vec(3*nz+1:4*nz)
        vz(:)    = eigen_vec(4*nz+1:bignz)
        call get_dens_from_pot(kcount)

        do i=1, nz
           write(20,fmt='(12(e22.15,x))') dble(drho(i)), dimag(drho(i)), &
                dble(dpres(i)), dimag(dpres(i)), &
                dble(dpot(i)), dimag(dpot(i)), &
                dble(vx(i)), dimag(vx(i)), &
                dble(vy(i)), dimag(vy(i)), &
                dble(vz(i)), dimag(vz(i))
        enddo
        endif         

     enddo

     if(maximize.eqv..true.) then
     loc    = maxloc(growth_fixb) 
     kcount = loc(1) 

      call fill_matrices(kcount,bcount)
      call eigenvalue_problem(kcount,bcount)

      write(10,fmt='(3(e22.15,x))') growth(kcount,bcount), freq(kcount,bcount), kvals(kcount)

        dpres(:) = eigen_vec(1:nz)
        dpot(:)  = eigen_vec(nz+1:2*nz)
        vx(:)    = eigen_vec(2*nz+1:3*nz)
        vy(:)    = eigen_vec(3*nz+1:4*nz)
        vz(:)    = eigen_vec(4*nz+1:bignz)
        call get_dens_from_pot(kcount)

        do i=1, nz
           write(20,fmt='(12(e22.15,x))') dble(drho(i)), dimag(drho(i)), &
                dble(dpres(i)), dimag(dpres(i)), &
                dble(dpot(i)), dimag(dpot(i)), &
                dble(vx(i)), dimag(vx(i)), &
                dble(vy(i)), dimag(vy(i)), &
                dble(vz(i)), dimag(vz(i))
        enddo
        endif

  enddo
  close(10)
  close(20)
  
end program viscsg

subroutine fill_matrices(kgrid, bgrid)
  use global
  implicit none
  integer :: i
  integer, intent(in) :: kgrid, bgrid
  real*8, parameter :: two_thirds = 2d0/3d0, one_third = 1d0/3d0, four_thirds = 4d0/3d0
  real*8 :: rho, dlogrho, dlogrhonu, k, bcool, alpha, dalpha, cs2, alpha_b, dlogrhonu_b, theta

  k     = kvals(kgrid)
  bcool = bvals(bgrid)
  Q3d   = Qarr(bgrid)
  theta = tirr 
 
  !fill the small matrices 
  do i=1, nz
     rho       = dens(i)
     dlogrho   = dlogdens(i)
     alpha     = alpha2d(i,bgrid)
     alpha_b   = bulkv*alpha 
     dalpha    = dalpha2d(i,bgrid)
     dlogrhonu = dlogdens(i) 
     if(alpha .gt. 0d0) dlogrhonu = dlogrhonu + dalpha/alpha
     dlogrhonu_b = dlogrhonu !because we assume bulk visc is const mult of shear visc 
     cs2       = csq(i) 

     !pressure/energy
     L11(i,:) = ( (gmma-1d0)*alpha*smallq**2d0*lambda*bgmma/cs2 - 1d0/bcool )*T(i,:)
     L12(i,:) = (gmma-1d0)*alpha*smallq**2d0*(1d0+mu)*(Q3d/rho)*( Tpp(i,:) - k*k*T(i,:) ) + (theta*cs2*Q3d/(bcool*bgmma*rho))*(Tpp(i,:) - k*k*T(i,:))
     L13(i,:) = -(ii*k*cs2*gmma/bgmma)*T(i,:)
     L14(i,:) = -2d0*(gmma-1d0)*ii*k*smallq*alpha*T(i,:)
     L15(i,:) = -cs2*(gmma*Tp_odd(i,:)/bgmma + dlogrho*T_odd(i,:) )
   
     R11(i,:) =  T(i,:)
     R12(i,:) =  0d0 
     R13(i,:) =  0d0 
     R14(i,:) =  0d0
     R15(i,:) =  0d0 

     !density/potential 
     L21(i,:) = 0d0
     L22(i,:) = 0d0
     L23(i,:) = -ii*k*T(i,:)
     L24(i,:) = 0d0 
     L25(i,:) = -dlogrho*T_odd(i,:) - Tp_odd(i,:)
   
     R21(i,:) =  0d0
     R22(i,:) =  (Q3d/rho)*( Tpp(i,:) - k*k*T(i,:) )
     R23(i,:) =  0d0 
     R24(i,:) =  0d0
     R25(i,:) =  0d0 

     !x momentum
     L31(i,:) = -ii*k*T(i,:)
     L32(i,:) = -ii*k*T(i,:)
     L33(i,:) = alpha*( Tpp(i,:) + dlogrhonu*Tp(i,:) - four_thirds*k*k*T(i,:) ) - alpha_b*k*k*T(i,:)
     L34(i,:) = 2d0*T(i,:)
     L35(i,:) = ii*k*alpha*( one_third*Tp_odd(i,:) + dlogrhonu*T_odd(i,:) ) + ii*k*alpha_b*Tp_odd(i,:)
   
     R31(i,:) =  0d0
     R32(i,:) =  0d0 
     R33(i,:) =  T(i,:)
     R34(i,:) =  0d0
     R35(i,:) =  0d0 
     
     !y mom eqn 
     L41(i,:) = -ii*k*smallq*alpha*lambda*bgmma*T(i,:)/cs2
     L42(i,:) = -ii*k*smallq*alpha*(1d0+mu)*(Q3d/rho)*( Tpp(i,:) - k*k*T(i,:) ) 
     L43(i,:) = (smallq-2d0)*T(i,:)
     L44(i,:) = alpha*( Tpp(i,:) + dlogrhonu*Tp(i,:) - k*k*T(i,:) )
     L45(i,:) = 0d0 
   
     R41(i,:) = 0d0 
     R42(i,:) = 0d0 
     R43(i,:) = 0d0  
     R44(i,:) = T(i,:)
     R45(i,:) = 0d0 
   
     !z mom eqn 
     L51(i,:) = -Tp(i,:) - dlogrho*T(i,:)
     L52(i,:) =  dlogrho*(Q3d*cs2/rho)*( Tpp(i,:) - k*k*T(i,:) ) - Tp(i,:) 
     L53(i,:) =  ii*k*alpha*(one_third*Tp(i,:) - two_thirds*dlogrhonu*T(i,:) ) + ii*alpha_b*k*( Tp(i,:) + dlogrhonu_b*T(i,:) )
     L54(i,:) =  0d0 
     L55(i,:) =  alpha*(four_thirds*Tpp_odd(i,:) + four_thirds*dlogrhonu*Tp_odd(i,:) - k*k*T_odd(i,:)) + alpha_b*( Tpp_odd(i,:) + dlogrhonu_b*Tp_odd(i,:) )
   
     R51(i,:) =  0d0 
     R52(i,:) =  0d0 
     R53(i,:) =  0d0 
     R54(i,:) =  0d0 
     R55(i,:) =  T_odd(i,:)
  enddo

  !boundary conditions at the top (for phi, vx, vy, vz)
  i = nz

  rho       = dens(i)
  dlogrho   = dlogdens(i)
  cs2       = csq(i) 

  if(vbc.eq.'wall')then
     L21(i,:) = 0d0 
     L22(i,:) = Tp(i,:) + k*T(i,:)
     L23(i,:) = 0d0
     L24(i,:) = 0d0
     L25(i,:) = 0d0
     
     R21(i,:) = 0d0 
     R22(i,:) = 0d0 
     R23(i,:) = 0d0 
     R24(i,:) = 0d0 
     R25(i,:) = 0d0 
  endif
  if(vbc.eq.'free')then
     L21(i,:) = 0d0 
     L22(i,:) = 0d0
     L23(i,:) = 0d0
     L24(i,:) = 0d0
     L25(i,:) = -rho*T_odd(i,:)/Q3d
     
     R21(i,:) = 0d0 
     R22(i,:) = Tp(i,:) + k*T(i,:) 
     R23(i,:) = 0d0 
     R24(i,:) = 0d0 
     R25(i,:) = 0d0 
  endif

  L31(i,:) = 0d0 
  L32(i,:) = 0d0 
  L33(i,:) = Tp(i,:)
  L34(i,:) = 0d0
  L35(i,:) = 0d0

  R31(i,:) = 0d0 
  R32(i,:) = 0d0 
  R33(i,:) = 0d0 
  R34(i,:) = 0d0 
  R35(i,:) = 0d0 

  L41(i,:) = 0d0 
  L42(i,:) = 0d0 
  L43(i,:) = 0d0 
  L44(i,:) = Tp(i,:)
  L45(i,:) = 0d0

  R41(i,:) = 0d0 
  R42(i,:) = 0d0 
  R43(i,:) = 0d0 
  R44(i,:) = 0d0 
  R45(i,:) = 0d0
  
  if(vbc.eq.'wall')then
     L51(i,:) = 0d0 
     L52(i,:) = 0d0 
     L53(i,:) = 0d0 
     L54(i,:) = 0d0
     L55(i,:) = T_odd(i,:)
     
     R51(i,:) = 0d0 
     R52(i,:) = 0d0 
     R53(i,:) = 0d0 
     R54(i,:) = 0d0 
     R55(i,:) = 0d0
  endif
  if(vbc.eq.'free')then
     L51(i,:) = 0d0 
     L52(i,:) = 0d0 
     L53(i,:) = 0d0 
     L54(i,:) = 0d0
     L55(i,:) =-cs2*dlogrho*T_odd(i,:)
     
     R51(i,:) = T(i,:)
     R52(i,:) = 0d0  
     R53(i,:) = 0d0 
     R54(i,:) = 0d0 
     R55(i,:) = 0d0
  endif
     
  

  !fill big matrix 
  !energy eqn
  bigLmatrix(1:nz,1:nz)         = L11 
  bigLmatrix(1:nz,nz+1:2*nz)    = L12 
  bigLmatrix(1:nz,2*nz+1:3*nz)  = L13
  bigLmatrix(1:nz,3*nz+1:4*nz)  = L14
  bigLmatrix(1:nz,4*nz+1:bignz) = L15
  
  bigRmatrix(1:nz,1:nz)         = R11 
  bigRmatrix(1:nz,nz+1:2*nz)    = R12 
  bigRmatrix(1:nz,2*nz+1:3*nz)  = R13
  bigRmatrix(1:nz,3*nz+1:4*nz)  = R14
  bigRmatrix(1:nz,4*nz+1:bignz) = R15

  !density/potential eqn 
  bigLmatrix(nz+1:2*nz,1:nz)         = L21 
  bigLmatrix(nz+1:2*nz,nz+1:2*nz)    = L22 
  bigLmatrix(nz+1:2*nz,2*nz+1:3*nz)  = L23
  bigLmatrix(nz+1:2*nz,3*nz+1:4*nz)  = L24
  bigLmatrix(nz+1:2*nz,4*nz+1:bignz) = L25
  
  bigRmatrix(nz+1:2*nz,1:nz)         = R21 
  bigRmatrix(nz+1:2*nz,nz+1:2*nz)    = R22 
  bigRmatrix(nz+1:2*nz,2*nz+1:3*nz)  = R23
  bigRmatrix(nz+1:2*nz,3*nz+1:4*nz)  = R24
  bigRmatrix(nz+1:2*nz,4*nz+1:bignz) = R25

  !vx eqn 
  bigLmatrix(2*nz+1:3*nz,1:nz)         = L31 
  bigLmatrix(2*nz+1:3*nz,nz+1:2*nz)    = L32 
  bigLmatrix(2*nz+1:3*nz,2*nz+1:3*nz)  = L33
  bigLmatrix(2*nz+1:3*nz,3*nz+1:4*nz)  = L34
  bigLmatrix(2*nz+1:3*nz,4*nz+1:bignz) = L35

  bigRmatrix(2*nz+1:3*nz,1:nz)         = R31 
  bigRmatrix(2*nz+1:3*nz,nz+1:2*nz)    = R32 
  bigRmatrix(2*nz+1:3*nz,2*nz+1:3*nz)  = R33
  bigRmatrix(2*nz+1:3*nz,3*nz+1:4*nz)  = R34
  bigRmatrix(2*nz+1:3*nz,4*nz+1:bignz) = R35

  !vy eqn 
  bigLmatrix(3*nz+1:4*nz,1:nz)         = L41 
  bigLmatrix(3*nz+1:4*nz,nz+1:2*nz)    = L42 
  bigLmatrix(3*nz+1:4*nz,2*nz+1:3*nz)  = L43
  bigLmatrix(3*nz+1:4*nz,3*nz+1:4*nz)  = L44
  bigLmatrix(3*nz+1:4*nz,4*nz+1:bignz) = L45

  bigRmatrix(3*nz+1:4*nz,1:nz)         = R41 
  bigRmatrix(3*nz+1:4*nz,nz+1:2*nz)    = R42 
  bigRmatrix(3*nz+1:4*nz,2*nz+1:3*nz)  = R43
  bigRmatrix(3*nz+1:4*nz,3*nz+1:4*nz)  = R44
  bigRmatrix(3*nz+1:4*nz,4*nz+1:bignz) = R45

  !vz eqn 
  bigLmatrix(4*nz+1:bignz,1:nz)         = L51 
  bigLmatrix(4*nz+1:bignz,nz+1:2*nz)    = L52 
  bigLmatrix(4*nz+1:bignz,2*nz+1:3*nz)  = L53
  bigLmatrix(4*nz+1:bignz,3*nz+1:4*nz)  = L54
  bigLmatrix(4*nz+1:bignz,4*nz+1:bignz) = L55
  
  bigRmatrix(4*nz+1:bignz,1:nz)         = R51 
  bigRmatrix(4*nz+1:bignz,nz+1:2*nz)    = R52 
  bigRmatrix(4*nz+1:bignz,2*nz+1:3*nz)  = R53
  bigRmatrix(4*nz+1:bignz,3*nz+1:4*nz)  = R54
  bigRmatrix(4*nz+1:bignz,4*nz+1:bignz) = R55

end subroutine fill_matrices

subroutine eigenvalue_problem(kgrid, bgrid)!fixed k and bcool
  use global
  implicit none 
  character*1, parameter :: jobvl = 'N', jobvr='V'
  integer, intent(in) :: kgrid, bgrid 
  integer :: i, info, loc(1)
  integer :: lwork 
  complex*16, allocatable :: work(:)
  complex*16 :: apha(bignz), bta(bignz), eigen(bignz)
  complex*16 :: vl(bignz,bignz), vr(bignz,bignz)
  real*8, parameter :: tol_im = 1d-9
  real*8 :: rwork(8*bignz), eigen_re, eigen_im
  
  lwork = 4*bignz
  allocate(work(lwork))

  !eigenvalue problem 
  call zggev (jobvl, jobvr, bignz, bigLmatrix, bignz, bigRmatrix, bignz, apha, bta, VL, & 
       bignz, VR, bignz, WORK, LWORK, RWORK, INFO) 
  if(info .ne. 0 ) print*, 'eigen failed?', info

  !compute eigenvalues and filter (get rid of decaying modes or modes growing too fast)
  do i=1, bignz
     eigen(i) = apha(i)/bta(i)
     eigen_re = dble(eigen(i))
     eigen_im = dimag(eigen(i))

     if((no_decay .eqv..true.).and.(dble(eigen_re).le.0d0)) eigen(i) = (0d0, 0d0)
     if((no_oscil .eqv..true.).and.(abs(eigen_im/eigen_re).gt.tol_im)) eigen(i) = (0d0, 0d0)
     
  enddo

  loc = maxloc(dble(eigen))
  growth(kgrid,bgrid) = dble(eigen(loc(1)))
  freq(kgrid,bgrid)   = dimag(eigen(loc(1)))
  eigen_vec           = vr(:,loc(1))

end subroutine eigenvalue_problem

subroutine get_dens_from_pot(kgrid)
  use global
  implicit none
  integer :: i
  integer, intent(in) :: kgrid
  real*8 :: k 
  complex*16 :: poisson(nz)
  
  k = kvals(kgrid)
  poisson = matmul(Tpp, dpot) - k*k*matmul(T, dpot)
  do i=1, nz 
     poisson(i) = poisson(i)*Q3d*csq(i)/dens(i)
  enddo
  drho = matmul(invT,poisson)

end subroutine get_dens_from_pot

subroutine chebyshev_poly(l, zbar, T_l, dT_l, d2T_l)
  implicit none
  real*8, intent(in) :: l, zbar
  real*8, intent(out):: T_l, dT_l, d2T_l
  real*8 :: t, lt, lsq
  lsq = l*l
  t = acos(zbar)
  lt = l*t
  T_l = cos(lt)
  if(abs(zbar).lt.1d0) then
     dT_l = l*sin(lt)/sin(t)
     d2T_l= -lsq*cos(lt)/sin(t)**2 + l*cos(t)*sin(lt)/sin(t)**3
  else
     dT_l =lsq
     d2t_l =lsq*(lsq-1d0)/3d0
     if(zbar.eq.-1d0) then
        dT_l = (-1d0)**(l+1d0)*dT_l
        d2T_l= (-1d0)**(l+2d0)*d2T_l
     endif
  endif
  return
end subroutine chebyshev_poly

subroutine get_den_and_dlogden(zout, d, dlogd)
  use global
  implicit none
  real*8, intent(in) :: zout
  real*8, intent(out):: d, dlogd
  integer, parameter :: neq = 2, itol = 1, itask=1, iopt=0, MF = 10 
  integer, parameter :: lrw =2*(20 + 16*NEQ), liw = 60
  integer            :: istate, iwork(liw), IPAR
  real*8, parameter  :: rtol = 1d-15, atol = 1d-15
  real*8             :: Y(neq), Z, rwork(lrw), RPAR, YDOT(neq)
  external           :: F, JAC
  
  Z    = 0d0 !start from midplane
  Y(1) = 1d0 !normalized density at z=0
  Y(2) = 0d0 !density gradient at z=0
  
  istate = 1 
  call DVODE (F, NEQ, Y, Z, ZOUT, ITOL, RTOL, ATOL, ITASK, &
       ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF, &
       RPAR, IPAR)
  !assign density and get log derivative
  d = Y(1)
  call F (NEQ, ZOUT, Y, YDOT, RPAR, IPAR)
  dlogd = YDOT(1)/d

end subroutine get_den_and_dlogden

subroutine get_thickness
  !iterate to get the height s.t. dens=0
  use global
  implicit none
  integer, parameter :: n = 1
  integer, parameter :: lwa = n*(3*n+13)
  integer            :: info 
  real*8,  parameter :: guess = 0.5d0, tol=1d-15
  real*8             :: x(n), fvec(n), wa(lwa)
  external           :: fcn 

  x(1) = guess

  call hybrd1(fcn,n,x,fvec,tol,info,wa,lwa)
  
  !physical height of domain boundary
  bigH = x(1)
  
end subroutine get_thickness

subroutine fcn(n,x,fvec,iflag)
  use global
  implicit none
  integer, parameter :: neq = 2, itol = 1, itask=1, iopt=0, MF = 10 
  integer, parameter :: lrw =2*(20 + 16*NEQ), liw = 60
  integer            :: n, iflag, istate, iwork(liw), IPAR
  real*8, parameter  :: rtol = 1d-15, atol = 1d-15
  real*8             :: x(n), fvec(n)
  real*8             :: Y(neq), Z, ZOUT, rwork(lrw), RPAR 
  real*8             :: final_density
  external           :: F, JAC
  
  Z    = 0d0 !start from midplane
  ZOUT = x(1)!guess end point 
  Y(1) = 1d0 !normalized density at z=0
  Y(2) = 0d0 !density gradient at z=0
  
  istate = 1 
  call DVODE (F, NEQ, Y, Z, ZOUT, ITOL, RTOL, ATOL, ITASK, &
       ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF, &
       RPAR, IPAR)

  !test if final density is small enough
  final_density = Y(1)
  fvec(1) = final_density - (1d0 - zmax)
  
end subroutine fcn

subroutine F (NEQ, Z, Y, YDOT, RPAR, IPAR)
  use global
  implicit none
  integer :: NEQ, IPAR
  real*8  :: Z, Y(NEQ), YDOT(NEQ), RPAR
  real*8  :: rho, gz

  rho = Y(1)
  gz  = Y(2)

  YDOT(1) = omegaz**2d0*Z + gz
  YDOT(1) = -YDOT(1)*rho**(2d0-bgmma)
  YDOT(2) = rho/Q3d

end subroutine F

SUBROUTINE JAC (NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
  DOUBLE PRECISION T, Y(NEQ), PD(NROWPD,NEQ), RPAR
end SUBROUTINE JAC

subroutine basic_setup(ncells)
  use global
  implicit none
  character*3 :: nbstring
  integer, intent(in) :: ncells
  integer :: i , j, jmid
  real*8 :: dz, zout, d, dlogd, amid, bcool
  real*8, external :: gviscosity

  !allocate arrays for basic state 
  allocate(zaxis(ncells))
  allocate(dens(ncells))
  allocate(dlogdens(ncells))
  allocate(csq(ncells))
  allocate(dcsq(ncells))
 
  !setup viscosity/heating array
  allocate(alpha2d(ncells,nb))
  allocate(dalpha2d(ncells,nb))
  allocate(Qarr(nb))

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
  
  !setup eqm density, sound-speed^2 and gradients 
  do i=1, ncells 
     zout = zaxis(i)
     call get_den_and_dlogden(zout, d, dlogd)
     dens(i)     = d
     dlogdens(i) = dlogd
     csq(i)      = d**(bgmma-1d0)
     
     dcsq(i)     = (bgmma-1d0)*csq(i)*dlogd
  enddo
  
  !setup viscosity/heating array 
  if(fixalpha .lt. 0d0) then !self-consistent heating/cooling
     do i=1, nb
        bcool = bvals(i)
        alpha2d(:,i)  = csq(:)/(bgmma*(gmma-1d0)*bcool*smallq**2d0)*(1d0 - tirr)
       dalpha2d(:,i) = dcsq(:)/(bgmma*(gmma-1d0)*bcool*smallq**2d0)*(1d0 - tirr)
     enddo
  else !viscosity, cooling and heating are independent (invoke unspecified heat source/sink) 
     do i=1, nb
        alpha2d(:,i)  = fixalpha
        dalpha2d(:,i) = 0d0
     enddo
  endif
 
  if(gvisc.eqv..true.) then !reset Q to depend on midplane alpha 
     do i=1, nb
     bcool = bvals(i) 
     amid    = 1d0/(bgmma*(gmma-1d0)*bcool*smallq**2d0)*(1d0 - tirr)
     Qarr(i) = gviscosity(amid) 
     enddo 
  else 
     Qarr(:) = Q3d 
  endif 
 
  if(ncells.eq.nzout) then !output basic state and deallocate arrays
     !output basic state 
     open(10,file='basic.dat')
     open(20,file='visc.dat')
     open(30,file='dvisc.dat')
     write(nbstring, fmt='(I3)') nb
     do i=1, nzout
        write(10,fmt='(5(e22.15,x))') zaxis(i), dens(i), dlogdens(i), csq(i), dcsq(i)
        write(20,fmt='('//trim(adjustl(nbstring))//'(e22.15,x))')  alpha2d(i,:)
        write(30,fmt='('//trim(adjustl(nbstring))//'(e22.15,x))') dalpha2d(i,:)
     enddo
     close(10)
     close(20)
     close(30)
     
     open(10,file='Qvals.dat')
     do i=1, nb
        write(10,fmt='(e22.15)') Qarr(i)
     enddo
     close(10)

     deallocate(zaxis)
     deallocate(dens)
     deallocate(dlogdens)
     deallocate(csq)
     deallocate(dcsq)
     deallocate(alpha2d)
     deallocate(dalpha2d)
     deallocate(Qarr)
     
  endif

end subroutine basic_setup

real*8 function gviscosity(avisc) 
use global
implicit none
real*8, intent(in) :: avisc 

gviscosity  = Q3d/sqrt(avisc) 
end function gviscosity 
