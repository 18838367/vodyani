program Q2
implicit none
real*8, allocatable, dimension(:) :: psi0, psi1, psi2, psi3, V, x
integer :: N0, N1, N2, N3, N, outunit, i
real*8 :: omega, dx, Emax, Emin, xmax, xmin, E

print*, "print this"
omega=1
dx=0.01
Emax=4
Emin=0
xmax=5
xmin=-5
N0=1
N1=2
N2=3
N3=4
E=3.5
N=ceiling((xmax-xmin)/dx)
allocate(psi0(N), psi1(N), psi2(N), psi3(N), V(N), x(N))
do i=1, N
    x(i)=dx*i+xmin
enddo

call potential(x,V,N,omega)
call numerov(V, psi0, dx, E, N)
open(newunit=outunit, file="Numerov.out", action="write")
do i=1, N
    write(outunit,*) psi0(i)
enddo
close(outunit)
end program Q2




subroutine numerov(V, psi, dx, E, N)
    implicit none
    integer, intent(in) :: N
    real*8, dimension(N), intent(in) :: V
    real*8, dimension(N), intent(out) :: psi
    real*8, intent(in) :: dx, E
    real*8, dimension(N) :: g
    integer :: i, m
    psi(1)=0
    psi(2)=0.0001
    psi(N)=0
    psi(N-1)=0.0001

    g=2*(V-E)
    m=ceiling(N/2.0)-50

    do i=2, m-1
        psi(i+1)=(2*(1+5*dx**2/12*g(i))*psi(i)-(1-dx**2/12*g(i-1))*psi(i-1))/(1-dx**2/12*g(i+1))
    enddo
    do i=2, N-m-1
        psi(N-i-1)=-1*(psi(N-i+1)*(1-dx**2/12*g(N-i+1))-2*(1+5*dx**2/12*g(N-i))*psi(N-i))/(1-dx**2/12*g(N-i-1))
    enddo
end subroutine numerov

subroutine potential(x,V,N,omega)
    implicit none
    integer, intent(in) :: N
    real*8, intent(in) :: omega
    real*8, dimension(N), intent(in) :: x
    real*8 , dimension(N), intent(out) :: V
    V=0.5*omega**2*x**2
end subroutine potential
