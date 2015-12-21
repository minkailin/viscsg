module global
  logical :: maximize, no_decay, no_oscil, betacool 
  character*4 :: vbc, tbc  
  integer, parameter :: nvar = 6
  integer :: nz, bignz, nr, nk, lmax,  nzout 
  real*8, parameter :: omegaz=1d0, pi = 2d0*acos(0d0)
  real*8, parameter :: stefan = 5.67d-5, bigG = 6.67d-8, mstar = 1.99d33, au=1.50d13, kboltz = 1.38d-16, mh=1.66d-24
  real*8, parameter :: mean_mu = 2.3, rstar = kboltz/mh, sig0=2200d0, T0=120d0, opacity_b2 = 5d-4, opacity_b0 = 0.24, year = 3.16d7
  real*8 :: fixQ, mdot, opscale, alpha_m, tirr, rmin, rmax, radius, alpha, opacity_law, opacity0, Tmid
  real*8, allocatable :: Q2d_arr(:), alpha_arr(:), raxis(:), temp_radial(:), bcool2d(:), theta2d(:), tau2d(:)

  real*8 :: Q3d, Q2d, smallq, gmma, bigH, mu, lambda, kmin, kmax, bulkv, zmax
  real*8, allocatable :: kvals(:) 
  real*8, allocatable :: zaxis(:), dens(:), dlogdens(:), pres(:), dlogpres(:), temp(:), dlogtemp(:), d2logtemp(:), sg(:), dsg(:), radflux(:)
  
  real*8, allocatable :: T(:,:), Tp(:,:), Tpp(:,:)
  real*8, allocatable :: T_odd(:,:), Tp_odd(:,:), Tpp_odd(:,:)
  real*8, allocatable :: growth(:,:), growth_fixr(:), freq(:,:)
  complex*16, parameter :: ii = (0d0, 1d0) 
  complex*16, allocatable :: drho(:), dpres(:), dpot(:), dtemp(:), vx(:), vy(:), vz(:) 
  complex*16, allocatable :: bigLmatrix(:,:), bigRmatrix(:,:), eigen_vec(:), eigen_fixr(:)
  complex*16, allocatable :: L11(:,:), L12(:,:), L13(:,:), L14(:,:), L15(:,:), L16(:,:)
  complex*16, allocatable :: L21(:,:), L22(:,:), L23(:,:), L24(:,:), L25(:,:), L26(:,:)
  complex*16, allocatable :: L31(:,:), L32(:,:), L33(:,:), L34(:,:), L35(:,:), L36(:,:)
  complex*16, allocatable :: L41(:,:), L42(:,:), L43(:,:), L44(:,:), L45(:,:), L46(:,:) 
  complex*16, allocatable :: L51(:,:), L52(:,:), L53(:,:), L54(:,:), L55(:,:), L56(:,:)
  complex*16, allocatable :: L61(:,:), L62(:,:), L63(:,:), L64(:,:), L65(:,:), L66(:,:)
  complex*16, allocatable :: R11(:,:), R12(:,:), R13(:,:), R14(:,:), R15(:,:), R16(:,:)
  complex*16, allocatable :: R21(:,:), R22(:,:), R23(:,:), R24(:,:), R25(:,:), R26(:,:)
  complex*16, allocatable :: R31(:,:), R32(:,:), R33(:,:), R34(:,:), R35(:,:), R36(:,:)
  complex*16, allocatable :: R41(:,:), R42(:,:), R43(:,:), R44(:,:), R45(:,:), R46(:,:)
  complex*16, allocatable :: R51(:,:), R52(:,:), R53(:,:), R54(:,:), R55(:,:), R56(:,:)
  complex*16, allocatable :: R61(:,:), R62(:,:), R63(:,:), R64(:,:), R65(:,:), R66(:,:)
end module global

program viscsg
  use global
  implicit none
  integer :: i,j, jmid, kcount, rcount, info, lwork
  integer :: loc(1) 
  integer, allocatable :: ipiv(:)
  real*8 :: dk, dlogr, optimize, dummy, eigen_re, eigen_im, kopt
  real*8 :: k, zbar, m, T_l, dT_l, d2T_l, n 
  real*8 :: surf, temperature, bcool, theta, thincool, omega 
  real*8, external :: toomreQ, bigOmega, optical_depth
  complex*16, allocatable :: work(:)
  complex*16 :: eigen
  namelist /basic2d/ fixQ, mdot, opscale, opacity_law, alpha_m, tirr, rmin, rmax, nr 
  namelist /params/ smallq, gmma, kmin, kmax, nk, mu, lambda, bulkv
  namelist /grid/   nz, zmax, nzout, vbc, tbc 
  namelist /job/    maximize, no_decay, no_oscil, betacool
  
  !read parameters
  open(7, file="input")
  read(7, nml=basic2d)
  read(7, nml=params)
  read(7, nml=grid)
  read(7, nml=job)
  close(7)

  if(opacity_law.eq.0d0) opacity0 = opacity_b0
  if(opacity_law.eq.2d0) opacity0 = opacity_b2

  mdot = mdot*mstar/year 

  !are we maximizing growth rates over k ?
  if(maximize.eqv..true.) then
     print*, 'maximize growth rates over k'
     optimize = 1d0
  else
    optimize= -1d0
 endif

 if(betacool.eqv..true.) then
    print*, 'setup vert iso disk, apply beta cool'
   thincool = 1d0
 else
    thincool = -1d0
 endif
 
  !output parameters for later analysis
  open(10,file='params.dat')
  write(10,fmt='(10(e22.15,x))') smallq, gmma, mu, lambda, tirr, dble(nz), optimize, bulkv, thincool, opacity_law
  close(10)
  
  !make sure grid cells are odd 
  if(mod(nz,2).eq.0) then
     print*, 'Nz needs to be odd but Nz=', nz
     stop
  endif
  
  !2D problem. setup radial grid and get alpha and Q2d at each radius 
  allocate(raxis(nr))
  allocate(Q2d_arr(nr))
  allocate(alpha_arr(nr))
  allocate(temp_radial(nr))
  allocate(bcool2d(nr))
  allocate(theta2d(nr))
  allocate(tau2d(nr))
  if(nr .gt. 1) then
     dlogr = log10(rmax/rmin)/(dble(nr)-1d0)
     do i=1, nr
        raxis(i) = 10d0**(log10(rmin) + dlogr*(dble(i)-1d0))
     enddo
  else
     raxis(1) = rmin 
  endif

  write(6,fmt='(A)') '************2D DISK************'
  
  surf        = 1d2
  temperature = 1d2
  do i=1, nr 
     radius = raxis(i)
     omega = bigOmega(radius)
     call get_surface_density_and_temperature(surf, temperature)
     call get_alpha(radius, surf, temperature)
     call get_bcool_and_theta(radius, surf, temperature, bcool, theta)
     alpha_arr(i)= alpha 
     Q2d_arr(i)  = toomreQ(radius, surf, temperature)
     bcool2d(i)     = bcool
     theta2d(i)     = theta 
!     if(betacool.eqv..true.) then !re-set temperature to that of optically thin cooling
!        temperature = alpha*smallq**2d0*omega*(rstar/mean_mu)
!        temperature = temperature/(4d0*stefan*opacity0*opscale)
!        temperature = temperature**(1d0/5d0)
!     endif
     temp_radial(i) = temperature 
     tau2d(i)       = optical_depth(surf, temperature)
  enddo

  !3D eigenvalue problem 

  allocate(kvals(nk))
  !wavenumber space
  if(nk .gt. 1) then
     dk = (kmax - kmin)/(dble(nk)-1d0)
     do i=1, nk 
        kvals(i) = kmin + dk*(dble(i)-1d0)
     enddo
  else
     kvals(1) = kmin
  endif
  open(10,file='kvals.dat')
  do i=1, nk
     write(10,fmt='(e22.15)') kvals(i)
  enddo
  close(10)
    
  allocate(T(nz,nz))
  allocate(Tp(nz,nz))
  allocate(Tpp(nz,nz))

  allocate(T_odd(nz,nz))
  allocate(Tp_odd(nz,nz))
  allocate(Tpp_odd(nz,nz))

!!$for the energy/pressure equation
  allocate(L11(nz,nz)) !operators on W
  allocate(L12(nz,nz)) !operators in temp 
  allocate(L13(nz,nz)) !operators on phi
  allocate(L14(nz,nz)) !operators on vx
  allocate(L15(nz,nz)) !operators on vy
  allocate(L16(nz,nz)) !operators on vz
  
  allocate(R11(nz,nz)) ; allocate(R12(nz,nz)) ; allocate(R13(nz,nz)) ; allocate(R14(nz,nz)) ; allocate(R15(nz,nz)) ; allocate(R16(nz,nz))

  !for the temperature/density equation (temperature definition plus continuity)
  allocate(L21(nz,nz)) ; allocate(L22(nz,nz)) ; allocate(L23(nz,nz)) ; allocate(L24(nz,nz)) ; allocate(L25(nz,nz)) ; allocate(L26(nz,nz))
  allocate(R21(nz,nz)) ; allocate(R22(nz,nz)) ; allocate(R23(nz,nz)) ; allocate(R24(nz,nz)) ; allocate(R25(nz,nz)) ; allocate(R26(nz,nz))

  !for density/potential equation (poisson plus continuity)
  allocate(L31(nz,nz)) ; allocate(L32(nz,nz)) ; allocate(L33(nz,nz)) ; allocate(L34(nz,nz)) ; allocate(L35(nz,nz)) ; allocate(L36(nz,nz))
  allocate(R31(nz,nz)) ; allocate(R32(nz,nz)) ; allocate(R33(nz,nz)) ; allocate(R34(nz,nz)) ; allocate(R35(nz,nz)) ; allocate(R36(nz,nz))

  !for the x momentum equation
  allocate(L41(nz,nz)) ; allocate(L42(nz,nz)) ; allocate(L43(nz,nz)) ; allocate(L44(nz,nz)) ; allocate(L45(nz,nz)) ; allocate(L46(nz,nz))
  allocate(R41(nz,nz)) ; allocate(R42(nz,nz)) ; allocate(R43(nz,nz)) ; allocate(R44(nz,nz)) ; allocate(R45(nz,nz)) ; allocate(R46(nz,nz))

  !for the y momentum equation
  allocate(L51(nz,nz)) ; allocate(L52(nz,nz)) ; allocate(L53(nz,nz)) ; allocate(L54(nz,nz)) ; allocate(L55(nz,nz)) ; allocate(L56(nz,nz))
  allocate(R51(nz,nz)) ; allocate(R52(nz,nz)) ; allocate(R53(nz,nz)) ; allocate(R54(nz,nz)) ; allocate(R55(nz,nz)) ; allocate(R56(nz,nz))

  !for the z  momentum equation
  allocate(L61(nz,nz)) ; allocate(L62(nz,nz)) ; allocate(L63(nz,nz)) ; allocate(L64(nz,nz)) ; allocate(L65(nz,nz)) ; allocate(L66(nz,nz))
  allocate(R61(nz,nz)) ; allocate(R62(nz,nz)) ; allocate(R63(nz,nz)) ; allocate(R64(nz,nz)) ; allocate(R65(nz,nz)) ; allocate(R66(nz,nz))

  bignz = nvar*nz
  allocate(bigLmatrix(bignz,bignz))
  allocate(bigRmatrix(bignz,bignz))
  allocate(growth(nk,nr))
  allocate(growth_fixr(nk))
  allocate(eigen_fixr(nk))
  allocate(freq(nk,nr))
  allocate(eigen_vec(bignz))
  allocate(drho(nz))
  allocate(dpres(nz))
  allocate(dpot(nz))
  allocate(dtemp(nz))
  allocate(vx(nz))
  allocate(vy(nz))
  allocate(vz(nz))


  open(10,file='basic2d.dat')
  open(20,file='basic3d.dat')
  open(30,file='eigenvalues.dat')
  open(40,file='eigenvectors.dat')
  open(50,file='output.dat',action='read') !2D data 

  write(6,fmt='(A)') '************3D DISK************'
  
  do rcount=1, nr 

     !need to solve for vertical structure at every radius 
     radius= raxis(rcount)
     Q2d   = Q2d_arr(rcount)
     Tmid  = temp_radial(rcount)
     
     if(betacool.eqv..true.) then
        alpha = 0d0 !we are going to use beta cooling in the pert eqns so set vert iso basic state 
     else
        alpha = alpha_arr(rcount) !self-consistent radiative cooling
     endif

     call get_thickness_and_Q3d
     
     write(6,fmt='(A,e22.15,e22.15,e22.15)') 'R, H, Q3d =', radius, bigH, Q3d
     
     write(10,fmt='(9(e22.15,x))') radius, bigH, Q3d, Q2d_arr(rcount), alpha_arr(rcount), bcool2d(rcount), theta2d(rcount), temp_radial(rcount), tau2d(rcount) 
     
     !setup basic state and output (fine grid for data analysis)
     call basic_setup(nzout)
     do j=1, nzout
        write(20,fmt='(11(e22.15,x))') zaxis(j), dens(j), dlogdens(j), pres(j), dlogpres(j), temp(j), dlogtemp(j), d2logtemp(j), sg(j), dsg(j), radflux(j) 
     enddo
     deallocate(zaxis)
     deallocate(dens)
     deallocate(dlogdens)
     deallocate(pres)
     deallocate(dlogpres)
     deallocate(temp)
     deallocate(dlogtemp)
     deallocate(d2logtemp)
     deallocate(sg)
     deallocate(dsg)
     deallocate(radflux)

     !setup basic state again (coarse spectral grid for eigen problem) 
     call basic_setup(nz) 
     !set up basis vectors 
     
     do j=1, nz !jth physical grid
        zbar = zaxis(j)/bigH
        do i=1, nz !ith basis
           
           !fill the even basis (these operate on W, temperature, potential, vx, vy)
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
     
     !eigen problem here 
     !loop through wavenumber space 

     do kcount=1, nk
        call fill_matrices(kcount, rcount)
        call eigenvalue_problem(kcount,rcount)
        
        growth_fixr(kcount) = growth(kcount,rcount) 
        eigen_fixr(kcount)  = dcmplx(growth(kcount,rcount), freq(kcount,rcount))         
 
        if(maximize.eqv..false.) then !we output each k 
           write(30,fmt='(3(e22.15,x))') growth(kcount,rcount), freq(kcount,rcount), kvals(kcount) 
           
           dpres(:) = eigen_vec(1:nz)
           dtemp(:) = eigen_vec(nz+1:2*nz)
           dpot(:)  = eigen_vec(2*nz+1:3*nz)
           vx(:)    = eigen_vec(3*nz+1:4*nz)
           vy(:)    = eigen_vec(4*nz+1:5*nz)
           vz(:)    = eigen_vec(5*nz+1:bignz)
           
           do i=1, nz
              write(40,fmt='(12(e22.15,x))') dble(dpres(i)), dimag(dpres(i)), &
                   dble(dtemp(i)), dimag(dtemp(i)), &
                   dble(dpot(i)), dimag(dpot(i)), &
                   dble(vx(i)), dimag(vx(i)), &
                   dble(vy(i)), dimag(vy(i)), &
                   dble(vz(i)), dimag(vz(i))
           enddo
        endif
     enddo
        
     if(maximize.eqv..true.) then !we output the most unstable k 
        loc    = maxloc(growth_fixr) 
!        read(50,fmt='(11(e22.15,x))') dummy, dummy, eigen_re, eigen_im, kopt, dummy, dummy, dummy, dummy, dummy, dummy 
!        eigen = dcmplx(eigen_re, eigen_im)
!        loc   = minloc(abs(eigen_fixr - eigen)**2d0 + abs(kvals-kopt)**2d0)

        kcount = loc(1)
        call fill_matrices(kcount, rcount)
        call eigenvalue_problem(kcount,rcount)
        !eigen_old = dcmplx(growth(kcount,rcount), freq(kcount,rcount))         

        write(30,fmt='(3(e22.15,x))') growth(kcount,rcount), freq(kcount,rcount), kvals(kcount)
        
        dpres(:) = eigen_vec(1:nz)
        dtemp(:) = eigen_vec(nz+1:2*nz)
        dpot(:)  = eigen_vec(2*nz+1:3*nz)
        vx(:)    = eigen_vec(3*nz+1:4*nz)
        vy(:)    = eigen_vec(4*nz+1:5*nz)
        vz(:)    = eigen_vec(5*nz+1:bignz)
           
        do i=1, nz
           write(40,fmt='(12(e22.15,x))') dble(dpres(i)), dimag(dpres(i)), &
                dble(dtemp(i)), dimag(dtemp(i)), &
                dble(dpot(i)), dimag(dpot(i)), &
                dble(vx(i)), dimag(vx(i)), &
                dble(vy(i)), dimag(vy(i)), &
                dble(vz(i)), dimag(vz(i))
        enddo
      
     endif
     
     
     !deallocate basic state after eigen calculation to prepare for next radius 
     deallocate(zaxis)
     deallocate(dens)
     deallocate(dlogdens)
     deallocate(pres)
     deallocate(dlogpres)
     deallocate(temp)
     deallocate(dlogtemp)
     deallocate(d2logtemp)     
     deallocate(sg)
     deallocate(dsg)
     deallocate(radflux)
  enddo

  close(10)
  close(20)
  close(30)
  close(40)  
  close(50)

end program viscsg

subroutine fill_matrices(kgrid, rcount)
  use global
  implicit none
  integer :: i
  integer, intent(in) :: kgrid, rcount
  real*8, parameter :: two_thirds = 2d0/3d0, one_third = 1d0/3d0, four_thirds = 4d0/3d0
  real*8 :: rho, dlogrho, dlogrhonu, k, alpha_b, dlogrhonu_b, fzero, rhomid, Omega, fzero_thin
  real*8 :: p, dlogp, tmp, dlogtmp, d2logtmp, bcool, theta 
  real*8, external :: bigOmega 

  k         = kvals(kgrid)
  alpha     = alpha_arr(rcount)
  bcool     = bcool2d(rcount)
  theta     = theta2d(rcount)
  Tmid      = temp_radial(rcount)

  alpha_b   = bulkv*alpha  !bulk visc is const mult of shear visc

  Omega  = bigOmega(radius)
  rhomid = Omega**2d0/(4d0*pi*bigG*Q3d)

  fzero  = 16d0*stefan*(gmma-1d0)*Omega
  fzero  = fzero/( 3d0*opacity0*opscale*rhomid**2d0*(rstar/mean_mu)**2d0 )
  fzero  = fzero*Tmid**(2d0 - opacity_law)

  fzero_thin = 4d0*stefan*(gmma-1d0)*opacity0*opscale*temp_radial(rcount)**6d0
  fzero_thin = fzero_thin/( (rstar/mean_mu)*temp_radial(rcount)*omega )

  !fill the small matrices 
  do i=1, nz
     rho       = dens(i)
     dlogrho   = dlogdens(i)
     p         = pres(i)
     dlogp     = dlogpres(i)
     tmp       = temp(i)
     dlogtmp   = dlogtemp(i)
     d2logtmp  = d2logtemp(i)
  
     dlogrhonu   = dlogrho   !because we assume alpha is constant wrt height 
     dlogrhonu_b = dlogrhonu !because we assume bulk visc is const mult of shear visc 
    
     !pressure/energy
     L11(i,:) = (gmma-1d0)*(alpha*smallq**2d0/tmp)*(1d0 + mu + lambda)*T(i,:)
     L12(i,:) =-(gmma-1d0)*(alpha*smallq**2d0/tmp)*(1d0 + mu)*T(i,:) 
     if(betacool.eqv..false.) then !radiative diffusion 
        L11(i,:) = L11(i,:) - fzero*(tmp**(3d0-opacity_law)/rho**2d0)*( ( (3d0-opacity_law)*dlogtmp - dlogrho )*dlogtmp*T(i,:) &
             + dlogtmp*Tp(i,:) + d2logtmp*T(i,:) )
        L12(i,:) = L12(i,:) - fzero*k*k*(tmp**(3d0-opacity_law)/rho**2d0)*T(i,:) &
             +fzero*(tmp**(3d0-opacity_law)/rho**2d0)*( ((3d0-opacity_law)*dlogtmp - dlogrho)*( Tp(i,:) + (4d0-opacity_law)*dlogtmp*T(i,:) ) &
             + Tpp(i,:) + (4d0-opacity_law)*dlogtmp*Tp(i,:) + (4d0-opacity_law)*d2logtmp*T(i,:) ) 
     else !beta cooling
        L11(i,:) = L11(i,:) - (1d0-theta)*T(i,:)/bcool
        L12(i,:) = L12(i,:) - theta*T(i,:)/bcool
        !proper optically thin coolling
        !        L11(i,:) =-T(i,:)
        !        L12(i,:) =-5d0*T(i,:)
     endif
       
     L13(i,:) = 0d0 
     L14(i,:) = -gmma*tmp*ii*k*T(i,:) 
     L15(i,:) = -(gmma-1d0)*2d0*alpha*smallq*ii*k*T(i,:)
     L16(i,:) = -tmp*dlogp*T_odd(i,:) - gmma*tmp*Tp_odd(i,:)

     R11(i,:) =  T(i,:)
     R12(i,:) =  0d0 
     R13(i,:) =  0d0 
     R14(i,:) =  0d0
     R15(i,:) =  0d0 
     R16(i,:) =  0d0 

     !density/temperature
     L21(i,:) = 0d0
     L22(i,:) = 0d0
     L23(i,:) = 0d0 
     L24(i,:) = -tmp*ii*k*T(i,:)
     L25(i,:) =  0d0 
     L26(i,:) = -tmp*( dlogrho*T_odd(i,:) + Tp_odd(i,:) )
   
     R21(i,:) =  T(i,:)
     R22(i,:) = -T(i,:)
     R23(i,:) =  0d0 
     R24(i,:) =  0d0
     R25(i,:) =  0d0 
     R26(i,:) =  0d0

     !density/potential 
     L31(i,:) = 0d0
     L32(i,:) = 0d0
     L33(i,:) = 0d0 
     L34(i,:) = -rho*ii*k*T(i,:)
     L35(i,:) =  0d0 
     L36(i,:) = -rho*( dlogrho*T_odd(i,:) + Tp_odd(i,:) )
     
     R31(i,:) =  0d0
     R32(i,:) =  0d0
     R33(i,:) =  Q3d*( Tpp(i,:) - k*k*T(i,:) )
     R34(i,:) =  0d0
     R35(i,:) =  0d0 
     R36(i,:) =  0d0
     
     !x momentum
     L41(i,:) = -ii*k*T(i,:)
     L42(i,:) = 0d0 
     L43(i,:) = -ii*k*T(i,:)
     L44(i,:) = alpha*( Tpp(i,:) + dlogrhonu*Tp(i,:) - four_thirds*k*k*T(i,:) ) - alpha_b*k*k*T(i,:)
     L45(i,:) = 2d0*T(i,:)
     L46(i,:) = ii*k*alpha*( one_third*Tp_odd(i,:) + dlogrhonu*T_odd(i,:) ) + ii*k*alpha_b*Tp_odd(i,:)
    
     R41(i,:) =  0d0
     R42(i,:) =  0d0 
     R43(i,:) =  0d0 
     R44(i,:) =  T(i,:)
     R45(i,:) =  0d0
     R46(i,:) =  0d0 
     
     !y mom eqn 
     L51(i,:) = -ii*k*smallq*(alpha/tmp)*(1d0 + mu + lambda)*T(i,:)
     L52(i,:) =  ii*k*smallq*(alpha/tmp)*(1d0 + mu)*T(i,:)
     L53(i,:) = 0d0 
     L54(i,:) = (smallq-2d0)*T(i,:)
     L55(i,:) = alpha*( Tpp(i,:) + dlogrhonu*Tp(i,:) - k*k*T(i,:) )
     L56(i,:) = 0d0 
   
     R51(i,:) = 0d0 
     R52(i,:) = 0d0 
     R53(i,:) = 0d0  
     R54(i,:) = 0d0 
     R55(i,:) = T(i,:)
     R56(i,:) = 0d0 
   
     !z mom eqn 
     L61(i,:) = -Tp(i,:) + (dlogp - dlogrho)*T(i,:)
     L62(i,:) = -dlogp*T(i,:)
     L63(i,:) = -Tp(i,:)
     L64(i,:) =  ii*k*alpha*(one_third*Tp(i,:) - two_thirds*dlogrhonu*T(i,:) ) + ii*alpha_b*k*( Tp(i,:) + dlogrhonu_b*T(i,:) )
     L65(i,:) =  0d0 
     L66(i,:) =  alpha*(four_thirds*Tpp_odd(i,:) + four_thirds*dlogrhonu*Tp_odd(i,:) - k*k*T_odd(i,:)) &
          + alpha_b*( Tpp_odd(i,:) + dlogrhonu_b*Tp_odd(i,:) )
   
     R61(i,:) =  0d0 
     R62(i,:) =  0d0 
     R63(i,:) =  0d0 
     R64(i,:) =  0d0 
     R65(i,:) =  0d0 
     R66(i,:) =  T_odd(i,:)
  enddo

  !boundary conditions at the top (for phi, vx, vy, vz)
  i = nz

  rho       = dens(i)
  dlogp     = dlogpres(i)
  tmp       = temp(i) 
  dlogtmp   = dlogtemp(i)

  !no temperature perturbation  (or no gradient?)
  if(betacool.eqv..false.) then !rad diffusion requires additional BC (here take no temp pert)
     L21(i,:) = 0d0 
     if(tbc.eq.'zero') L22(i,:) = T(i,:) 
     if(tbc.eq.'flat') L22(i,:) = Tp(i,:)
     
     !no pert z flux
!     L21(i,:)  =-dlogtmp*T(i,:)
!     L22(i,:)  = Tp(i,:) + 2d0*dlogtmp*T(i,:)

     L23(i,:) = 0d0
     L24(i,:) = 0d0
     L25(i,:) = 0d0
     L26(i,:) = 0d0
     
     R21(i,:) = 0d0 
     R22(i,:) = 0d0 
     R23(i,:) = 0d0 
     R24(i,:) = 0d0 
     R25(i,:) = 0d0 
     R26(i,:) = 0d0 
  endif

  !potential condition 
  if(vbc.eq.'wall')then
     L31(i,:) = 0d0 
     L32(i,:) = 0d0 
     L33(i,:) = Tp(i,:) + k*T(i,:)
     L34(i,:) = 0d0
     L35(i,:) = 0d0
     L36(i,:) = 0d0
     
     R31(i,:) = 0d0 
     R32(i,:) = 0d0 
     R33(i,:) = 0d0 
     R34(i,:) = 0d0 
     R35(i,:) = 0d0 
     R36(i,:) = 0d0
  endif
  if(vbc.eq.'free')then
     L31(i,:) = 0d0 
     L32(i,:) = 0d0
     L33(i,:) = 0d0
     L34(i,:) = 0d0
     L35(i,:) = 0d0 
     L36(i,:) = -rho*T_odd(i,:)/Q3d
     
     R31(i,:) = 0d0 
     R32(i,:) = 0d0 
     R33(i,:) = Tp(i,:) + k*T(i,:) 
     R34(i,:) = 0d0 
     R35(i,:) = 0d0 
     R36(i,:) = 0d0 
  endif

  !no vx grad 
  L41(i,:) = 0d0 
  L42(i,:) = 0d0 
  L43(i,:) = 0d0
  L44(i,:) = Tp(i,:)
  L45(i,:) = 0d0
  L46(i,:) = 0d0 

  R41(i,:) = 0d0 
  R42(i,:) = 0d0 
  R43(i,:) = 0d0 
  R44(i,:) = 0d0 
  R45(i,:) = 0d0 
  R46(i,:) = 0d0 

  !no vy grad 
  L51(i,:) = 0d0 
  L52(i,:) = 0d0 
  L53(i,:) = 0d0 
  L54(i,:) = 0d0
  L55(i,:) = Tp(i,:)
  L56(i,:) = 0d0 

  R51(i,:) = 0d0 
  R52(i,:) = 0d0 
  R53(i,:) = 0d0 
  R54(i,:) = 0d0 
  R55(i,:) = 0d0
  R56(i,:) = 0d0 
  
  if(vbc.eq.'wall')then
     L61(i,:) = 0d0 
     L62(i,:) = 0d0 
     L63(i,:) = 0d0 
     L64(i,:) = 0d0
     L65(i,:) = 0d0
     L66(i,:) = T_odd(i,:)
      
     R61(i,:) = 0d0 
     R62(i,:) = 0d0 
     R63(i,:) = 0d0 
     R64(i,:) = 0d0 
     R65(i,:) = 0d0
     R66(i,:) = 0d0
  endif
  if(vbc.eq.'free')then
     L61(i,:) = 0d0 
     L62(i,:) = 0d0 
     L63(i,:) = 0d0 
     L64(i,:) = 0d0
     L65(i,:) = 0d0
     L66(i,:) =-tmp*dlogp*T_odd(i,:)
     
     R61(i,:) = T(i,:)
     R62(i,:) = 0d0  
     R63(i,:) = 0d0 
     R64(i,:) = 0d0 
     R65(i,:) = 0d0
     R66(i,:) = 0d0
  endif
     
  !fill big matrix 

  !energy eqn
  bigLmatrix(1:nz,1:nz)         = L11 
  bigLmatrix(1:nz,nz+1:2*nz)    = L12 
  bigLmatrix(1:nz,2*nz+1:3*nz)  = L13
  bigLmatrix(1:nz,3*nz+1:4*nz)  = L14
  bigLmatrix(1:nz,4*nz+1:5*nz)  = L15
  bigLmatrix(1:nz,5*nz+1:bignz) = L16

  bigRmatrix(1:nz,1:nz)         = R11 
  bigRmatrix(1:nz,nz+1:2*nz)    = R12 
  bigRmatrix(1:nz,2*nz+1:3*nz)  = R13
  bigRmatrix(1:nz,3*nz+1:4*nz)  = R14
  bigRmatrix(1:nz,4*nz+1:5*nz)  = R15
  bigRmatrix(1:nz,5*nz+1:bignz) = R16

  !temperature/density 
  bigLmatrix(nz+1:2*nz,1:nz)         = L21 
  bigLmatrix(nz+1:2*nz,nz+1:2*nz)    = L22 
  bigLmatrix(nz+1:2*nz,2*nz+1:3*nz)  = L23
  bigLmatrix(nz+1:2*nz,3*nz+1:4*nz)  = L24
  bigLmatrix(nz+1:2*nz,4*nz+1:5*nz)  = L25
  bigLmatrix(nz+1:2*nz,5*nz+1:bignz) = L26

  bigRmatrix(nz+1:2*nz,1:nz)         = R21 
  bigRmatrix(nz+1:2*nz,nz+1:2*nz)    = R22 
  bigRmatrix(nz+1:2*nz,2*nz+1:3*nz)  = R23
  bigRmatrix(nz+1:2*nz,3*nz+1:4*nz)  = R24
  bigRmatrix(nz+1:2*nz,4*nz+1:5*nz)  = R25
  bigRmatrix(nz+1:2*nz,5*nz+1:bignz) = R26

  !density/potential 
  bigLmatrix(2*nz+1:3*nz,1:nz)         = L31 
  bigLmatrix(2*nz+1:3*nz,nz+1:2*nz)    = L32 
  bigLmatrix(2*nz+1:3*nz,2*nz+1:3*nz)  = L33
  bigLmatrix(2*nz+1:3*nz,3*nz+1:4*nz)  = L34
  bigLmatrix(2*nz+1:3*nz,4*nz+1:5*nz)  = L35
  bigLmatrix(2*nz+1:3*nz,5*nz+1:bignz) = L36

  bigRmatrix(2*nz+1:3*nz,1:nz)         = R31 
  bigRmatrix(2*nz+1:3*nz,nz+1:2*nz)    = R32 
  bigRmatrix(2*nz+1:3*nz,2*nz+1:3*nz)  = R33
  bigRmatrix(2*nz+1:3*nz,3*nz+1:4*nz)  = R34
  bigRmatrix(2*nz+1:3*nz,4*nz+1:5*nz)  = R35
  bigRmatrix(2*nz+1:3*nz,5*nz+1:bignz) = R36

  !vx eqn 
  bigLmatrix(3*nz+1:4*nz,1:nz)         = L41 
  bigLmatrix(3*nz+1:4*nz,nz+1:2*nz)    = L42 
  bigLmatrix(3*nz+1:4*nz,2*nz+1:3*nz)  = L43
  bigLmatrix(3*nz+1:4*nz,3*nz+1:4*nz)  = L44
  bigLmatrix(3*nz+1:4*nz,4*nz+1:5*nz)  = L45
  bigLmatrix(3*nz+1:4*nz,5*nz+1:bignz) = L46

  bigRmatrix(3*nz+1:4*nz,1:nz)         = R41 
  bigRmatrix(3*nz+1:4*nz,nz+1:2*nz)    = R42 
  bigRmatrix(3*nz+1:4*nz,2*nz+1:3*nz)  = R43
  bigRmatrix(3*nz+1:4*nz,3*nz+1:4*nz)  = R44
  bigRmatrix(3*nz+1:4*nz,4*nz+1:5*nz)  = R45
  bigRmatrix(3*nz+1:4*nz,5*nz+1:bignz) = R46

  !vy eqn 
  bigLmatrix(4*nz+1:5*nz,1:nz)         = L51 
  bigLmatrix(4*nz+1:5*nz,nz+1:2*nz)    = L52 
  bigLmatrix(4*nz+1:5*nz,2*nz+1:3*nz)  = L53
  bigLmatrix(4*nz+1:5*nz,3*nz+1:4*nz)  = L54
  bigLmatrix(4*nz+1:5*nz,4*nz+1:5*nz)  = L55
  bigLmatrix(4*nz+1:5*nz,5*nz+1:bignz) = L56

  bigRmatrix(4*nz+1:5*nz,1:nz)         = R51 
  bigRmatrix(4*nz+1:5*nz,nz+1:2*nz)    = R52 
  bigRmatrix(4*nz+1:5*nz,2*nz+1:3*nz)  = R53
  bigRmatrix(4*nz+1:5*nz,3*nz+1:4*nz)  = R54
  bigRmatrix(4*nz+1:5*nz,4*nz+1:5*nz)  = R55
  bigRmatrix(4*nz+1:5*nz,5*nz+1:bignz) = R56

  !vz eqn 
  bigLmatrix(5*nz+1:bignz,1:nz)         = L61 
  bigLmatrix(5*nz+1:bignz,nz+1:2*nz)    = L62 
  bigLmatrix(5*nz+1:bignz,2*nz+1:3*nz)  = L63
  bigLmatrix(5*nz+1:bignz,3*nz+1:4*nz)  = L64
  bigLmatrix(5*nz+1:bignz,4*nz+1:5*nz)  = L65
  bigLmatrix(5*nz+1:bignz,5*nz+1:bignz) = L66
  
  bigRmatrix(5*nz+1:bignz,1:nz)         = R61 
  bigRmatrix(5*nz+1:bignz,nz+1:2*nz)    = R62 
  bigRmatrix(5*nz+1:bignz,2*nz+1:3*nz)  = R63
  bigRmatrix(5*nz+1:bignz,3*nz+1:4*nz)  = R64
  bigRmatrix(5*nz+1:bignz,4*nz+1:5*nz)  = R65
  bigRmatrix(5*nz+1:bignz,5*nz+1:bignz) = R66

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
  real*8, parameter :: tol_im = 1d-9, maxrate=1d1
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

     if(abs(eigen(i)) .gt. maxrate ) eigen(i) = (0d0, 0d0)
     if((no_decay .eqv..true.).and.(eigen_re.le.0d0)) eigen(i) = (0d0, 0d0)
     if((no_oscil .eqv..true.).and.(abs(eigen_im/eigen_re).gt.tol_im)) eigen(i) = (0d0, 0d0)
     
  enddo

  loc = maxloc(dble(eigen))
  growth(kgrid,bgrid) = dble(eigen(loc(1)))
  freq(kgrid,bgrid)   = dimag(eigen(loc(1)))
  eigen_vec           = vr(:,loc(1))

end subroutine eigenvalue_problem


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

