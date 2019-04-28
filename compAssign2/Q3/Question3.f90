program Q3
implicit none
real*8 :: dR, k, rmax, pi
integer :: l, N, i, lambda
real*8, allocatable, dimension(:) :: rb, num, r
real*8, external :: riccati_bessel
rmax=10
dR=0.0001
pi=3.14159
N=rmax/dR
l=2
k=4

allocate(r(N), rb(N), num(N))
do i=1, N
    r(i)=i*dR
    rb(i)=riccati_bessel(l,r(i)*k)
enddo

call numerovForw(num, dR, N, r, l, k)
lambda=ceiling(2*pi/k/dR)
num=num/(maxval(num((N-lambda):N)))

open(unit=100, file="Q3.out", action="write")
do i=1, N
    write(100,*) r(i), rb(i), num(i)
enddo
close(100)
print*, "done"
end program Q3


subroutine numerovForw(psi, dx, N, r, l, k)
    implicit none
    integer, intent(in) :: N, l
    real*8, dimension(N), intent(in) :: r
    real*8, dimension(N), intent(inout) :: psi
    real*8, intent(in) :: dx, k
    real*8, dimension(N) :: g
    integer :: i
    psi(1)=0
    psi(2)=0.0001
    g=(l*(l+1))/(r**2)-k**2
    
    do i=2, N-1
        psi(i+1)=(2*(1+5*dx**2/12*g(i))*psi(i)-(1-dx**2/12*g(i-1))*psi(i-1))/(1-dx**2/12*g(i+1))
    enddo
end subroutine numerovForw

real*8 function riccati_bessel(l,x) result(S)
implicit none
integer, intent(in) :: l
real*8, intent(in) :: x
integer, external :: doubleFac
real*8 :: S0, S1
integer :: j

if(x<1) then
    S=(x**(l+1))/doubleFac(2*l+1)
else
    S0 = sin(x); S1 = sin(x)/x - cos(x)
    if(l==0) then
    S = S0
    elseif(l==1) then
    S = S1
    else
    do j=2,l
    S = dble(2*j-1)/x * S1 - S0
    S0 = S1; S1 = S
    enddo
    endif
endif
end function riccati_bessel


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
