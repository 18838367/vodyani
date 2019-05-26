!  This file is an outline of a code to perform a calculation of a
!    charged particle scattered by a spherically-symmetric short-range local
!    potential.
!  There are a number of comments throughout the file beginning with !>>> which
!    indicate sections of code you will need to complete.

program main

  use constants
  implicit none

  real*8, allocatable :: &
    Vmat(:,:),       & !V-matrix elements Vmat(kf,ki)
    V(:),            & !radial potential function V(r)
    rgrid(:),        & !radial grid
    rweights(:),     & !radial integration weights
    kgrid(:),        & !momentum-space grid
    kweights(:),     & !momentum-space integration weights (and Green's func)
    contwaves(:,:),  & !projectile radial continuum waves contwaves(r,k)
    DCS(:),          & !array to hold differential cross section - DCS(theta)
    theta(:),        & !array to hold values of theta - in degrees
    ICS(:)             !integrated cross section per l

  real*8 :: &
    rmax,   & !max value of radial grid
    dr,     & !radial grid step size
    energy, & !projectile energy
    k,      & !projectile momentum
    kg_A, kg_B, kg_P !some parameters for setting up kgrid

  complex*16, allocatable :: Ton(:) !on-shell T-matrix element per l

  integer :: &
    i,          & !do loop integer
    nrmax,      & !number of rgrid points
    nkmax,      & !number of kgrid points
    zproj,      & !projectile charge 
    l,          & !partial-wave angular momentum
    lmin, lmax, & !min and max values of l
    iounit,     & !a unit number for input/ouput
    ntheta,     & !an index to iterate over theta
    nthetamax,  & !max number of theta
    kg_Na, kg_Nb, kg_Np !some more parameters for setting up kgrid
  !set kgrid parameters 
  ! leave this as is, or modify to read them in from input file
  ! if you want to play around with it and understand how the
  ! kgrids are set up
    kg_Na = 30; kg_Nb = 30; kg_Np = 10
    kg_a = 0.85; kg_b = 2.5; kg_p = 4.0
    nkmax=kg_Na+kg_Nb+kg_Np+1

  !>>> open data.in file and read input parameters
  !    note: energy should be read in electron volts
  !      and grid parameters in atomic units
    open(unit=100, file="data.in", action="read")
    read(100,*) energy
    read(100,*) rmax, dr
    read(100,*) zproj, lmin, lmax
    close(100)

  !>>> do any input validation you think is necessary here
 
  !>>> convert the energy to atomic units and calculate the
  !      projectile momentum
    energy=energy*0.0367493
    k=(2.0d0*energy)**(0.5)
  !>>> determine number of rgrid points nrmax
  !    note: nrmax should be odd for simpson's integration
    nrmax=floor(rmax/dr)+mod(rmax,dr)
  !allocate memory
    allocate(rgrid(nrmax),rweights(nrmax))
    allocate(kgrid(nkmax),kweights(nkmax))
    allocate(contwaves(nrmax,nkmax))
    allocate(V(nrmax))
    allocate(Ton(lmin:lmax),ICS(lmin:lmax))
    allocate(Vmat(nkmax,nkmax))

  !setup grids
    call setup_rgrid(nrmax, dr, rgrid, rweights)
    call setup_kgrid(k, nkmax, kg_Na, kg_a, kg_Nb, kg_b, kg_Np, kg_p, kgrid, kweights)
  
  !>>> define short-range potential V(r)
    V=zproj*(1+1.0d0/(rgrid))*exp(-2.0d0*rgrid)
  !begin loop over angular momenta
  do l=lmin, lmax
    !populate contwaves matrix with a continuum wave for each off-shell k
      call setup_contwaves(nkmax,kgrid,l,nrmax,rgrid,contwaves)
    
    !evaluate the V-matrix elements  
      call calculate_Vmatrix(nkmax,kgrid,contwaves,nrmax,rgrid,rweights,V,Vmat)

    !solve the Lippman-Schwinger equation for the on-shell T-matrix
      call tmatrix_solver(nkmax,kgrid,kweights,Vmat,Ton(l))
  enddo

  !populate theta - theta = 0, 1,..., 180 degrees
  nthetamax = 181
  allocate(theta(nthetamax),DCS(nthetamax))
  do ntheta = 1, nthetamax
    theta(ntheta) = dble(ntheta-1)
  enddo

  !call subroutines to calculate DCS, ICS, and phase shifts
    call compute_dcs(nthetamax, theta, lmin, lmax, Ton, k, DCS)
    call compute_ics(lmin, lmax, Ton, k, ICS)

  !>>> output the DCS to file as a function of theta
    open(unit=200, file="DCS.out", action="write")
    do i=1, nrmax, 1
        write(200,*) theta(i), DCS(i)
    enddo
    close(200)
  !>>> either write the ICS to a file or print to screen
  !    in a way you can easily extract it using GREP
    open(unit=300, file="ICS.out", action="write")
    do i=1, nrmax, 1
        write(300,*) ICS(i) 
    enddo
    close(300)
end program

subroutine compute_ics(lmin, lmax, Ton, k, ICS)
  use constants
  implicit none
  integer, intent(in) :: lmin, lmax
  complex*16, intent(in) :: Ton(lmin:lmax)
  real*8, intent(in) :: k
  real*8, intent(out) :: ICS(lmin:lmax)
  integer :: l

  !>>> populate the ICS array with the partial-wave ICS per l
    ICS=4*3.14159**3/(k**4)*(2*l+1)*(abs(Ton)**2)
end subroutine compute_ics

subroutine compute_dcs(nthetamax, theta, lmin, lmax, Ton, k, DCS)
  use constants
  implicit none
  integer, intent(in) :: nthetamax, lmin, lmax
  real*8, intent(in) :: theta(nthetamax), k
  complex*16, intent(in) :: Ton(0:lmax)
  real*8, intent(out) :: DCS(nthetamax)
  integer :: l, ntheta !loop indices
  real*8 :: PL !Legendre polynomials - from file plql.f
  real*8 :: costheta !use this to store cos(theta in radians)
  complex*16 :: f(nthetamax) !scattering amplitude

  !>>> calculate the scattering amplitude f(theta) for each theta
  !    by iterating over l and using the partial-wave
  !    expansion of f
    f=0
    
    do ntheta=1, nthetamax, 1
        print*, "dont kill me..."
        print*, ntheta
        do l=lmin, lmax, 1
            print*, -3.14159/(k**2)*(2*l+1*Ton(l))*PL(costheta)
            costheta=cos(theta(ntheta))
            f(ntheta)=f(ntheta)-3.14159/(k**2)*(2*l+1*Ton(l))*PL(costheta)
        enddo
    enddo
    print*, "lol wut"
  !>>> obtain the DCS from the scattering amplitude
    DCS=abs(f)**2
end subroutine compute_dcs

subroutine setup_rgrid(nrmax, dr, rgrid, rweights)
  implicit none
  integer, intent(in) :: nrmax
  real*8, intent(in) :: dr
  real*8, intent(out) :: rgrid(nrmax), rweights(nrmax)
  integer :: ir !index to iterate over r

  !>>> iterate over r and populate the rgrid and rweights arrays
  !      - rweights should contain Simpson's integration weights:
  !        (4, 2, 4, 2, ..., 2, 4, 1) * dr / 3.0
  !      - you can make use of the intrinsic MOD function for the 
  !        alternating 4, 2 terms
  !      - note the first weight of 1 has been skipped as we skip the point r=0

    do ir=1, nrmax, 1
        rgrid(ir)=dr*ir
        rweights(ir)=2+2*MOD(ir, 2)
    enddo
    rweights(nrmax)=1
    rweights=rweights*dr/3.0

end subroutine setup_rgrid

subroutine setup_contwaves(nkmax, kgrid, l, nrmax, rgrid, contwaves)
  implicit none
  integer, intent(in) :: nkmax, l, nrmax
  real*8, intent(in) :: kgrid(nkmax), rgrid(nrmax)
  real*8, intent(out) :: contwaves(nrmax,nkmax)
  integer :: nk, nr !indices to loop over k and r
  real*8, external :: riccati_bessel !external function (provided in this file)
                                     !to evaluate the Riccati-Bessel functions
  
  !>>> iterate over k and r, populating the contwaves matrix
  !    (you may wish to parallelise this operation)
    do nk=1, nkmax, 1
        do nr=1, nrmax, 1 
            contwaves(nr,nk)=riccati_bessel(l,kgrid(nk)*rgrid(nr))
        enddo
    enddo
end subroutine setup_contwaves
    
subroutine calculate_Vmatrix(nkmax,kgrid,contwaves,nrmax,rgrid,rweights,V,Vmat)
  use constants
  implicit none
  integer, intent(in) :: nkmax, nrmax
  real*8, intent(in) :: kgrid(nkmax), contwaves(nrmax,nkmax), rgrid(nrmax), rweights(nrmax), V(nrmax)
  real*8, intent(out) :: Vmat(nkmax,nkmax)
  integer :: nkf,nki !indices for looping over on- and off-shell k

  !>>> evaluate the V-matrix elements and store in the Vmat matrix
  !    note: the V-matrix is symmetric, make use of this fact to reduce the 
  !          amount of time spent in this subroutine
  !    (you may wish to parallelise this operation)
    do nkf=1, nkmax, 1
        do nki=1, nkf, 1
            Vmat(nkf, nki)=sum(contwaves(:,nkf)*V*contwaves(:,nki)*rweights)*2/3.14159
            Vmat(nki, nkf)=Vmat(nkf, nki)        
        enddo
    enddo
    
end subroutine calculate_Vmatrix
    
subroutine tmatrix_solver(nkmax,kgrid,kweights,Vmat,Ton)
  use constants
  implicit none
  integer, intent(in) :: nkmax
  real*8, intent(in) :: kgrid(nkmax), kweights(nkmax), Vmat(nkmax,nkmax)
  complex*16, intent(out) :: Ton !on-shell T-matrix element
  real*8 :: &
    Koff(nkmax-1), & !half-off-shell K-matrix elements
    Kon,           & !on-shell K-matrix element
    Von,           & !on-shell V-matrix element
    A(nkmax-1,nkmax-1) !Coefficient matrix for the linear system Ax=b
  integer, external :: del
  integer :: i, j, ipiv(nkmax-1), info

  !>>> set up the coefficient matrix A for the linear system
    do i=1,nkmax-1, 1
        do j=1,nkmax-1, 1
            A(i,j)=del(i,j)-kweights(j+1)*Vmat(i+1,j+1)
        enddo
    enddo    
    
 
  !>>> store the RHS of the linear system in Koff
    do j=1, nkmax-1, 1
      Koff(j)=Vmat(j+1,1)
    enddo
  call dgesv( nkmax-1, 1, A, nkmax-1, ipiv, Koff, nkmax-1, info )
  if(info /= 0) then
    print*, 'ERROR in dgesv: info = ', info
  endif

  !>>> calculate the on-shell K-matrix element (Kon)
    Kon=Koff(1)
  !>>> calculate the on-shell T-matrix element (Ton)
    Ton=Kon/(1+CMPLX(0,1)*3.14159*Kon/kgrid(nkmax))
end subroutine tmatrix_solver

real*8 function riccati_bessel(l,x) result(RB)
implicit none
  !>>> Insert your code here for calculating the Riccati-Bessel functions
integer, intent(in) :: l
real*8, intent(in) :: x
integer, external :: doubleFac
real*8 :: S0, S1
integer :: j

if(x<1) then !if small argument use limit instead, checked with some diagnostic graphs for limit
    RB=(x**(l+1))/doubleFac(2*l+1)
else
    S0 = sin(x); S1 = sin(x)/x - cos(x)
    if(l==0) then
    RB = S0
    elseif(l==1) then
    RB = S1
    else
    do j=2,l
    RB = dble(2*j-1)/x * S1 - S0
    S0 = S1; S1 = RB
    enddo
    endif
endif
end function riccati_bessel

!calculates the double factial for input integer n
integer function doubleFac(n) result(f)
implicit none
integer, intent(in) :: n
integer :: i
i=n
f=1
do while(i>1)
    f=f*i
    i=i-2
enddo
end function doubleFac

!provides behaviour of a delta matrix (identity matrix)
integer function del(i,j) result(f)
implicit none
integer, intent(in) :: i, j

if(i==j) then !i.e. if not on the diagonal then zero
    f=1
else 
    f=0
endif
end function 
   

!A subroutine provided for you to set up the kgrid and kweights
!note: the kgrid is setup with the on-shell point in the first element
!      and the corresponding kweights include the integration weights
!      AND the Green's function
subroutine setup_kgrid(k,nkmax,Na,a,Nb,b,Np,p,kgrid,kweights)
  implicit none
  real*8, intent(in) :: k
  integer, intent(in) :: Na, Nb, Np
  integer, intent(in) :: nkmax
  real*8, intent(out) :: kgrid(nkmax), kweights(nkmax)
  integer :: nk
  real*8, intent(in) ::  a, b, p
  real*8 :: grid1(nkmax-1), weight1(nkmax-1)

  call kgrid_igor(0.0,k,a,Na,b,Nb,p,Np,nkmax-1,grid1,weight1)
  
  kgrid(1) = k
  kgrid(2:nkmax) = grid1
  kweights(1) = 0.0d0
  kweights(2:nkmax) = weight1
  do nk=2, nkmax
    kweights(nk) = 2.0d0 * kweights(nk) / (k**2 - kgrid(nk)**2)
  enddo
end subroutine setup_kgrid
