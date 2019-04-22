program bisection
implicit none
real*8 :: psi1, psi2, E, Emax, Emin, dx, xmax, xmin, omega, tol
real*8, allocatable, dimension(:) :: x, psi, V
integer :: i, N, outunit
dx=0.1
omega=1
Emin=0
Emax=4
xmin=-5
xmax=5
psi1=0
psi2=0.0001
N=ceiling((xmax-xmin)/dx)
allocate(x(N), psi(N), V(N))
do i=1, N
    x(i)=dx*(i-1)+xmin
enddo
call potential(x,V,N,omega)
psi(1)=psi1
psi(2)=psi2
tol=0.001


print*, E
E=(Emin+Emax)/2
call shoot(x, V, N, E, dx, psi)
if (psi(N)>0) then
    Emin=E
elseif (psi(N)<0) then 
    Emax=E
endif



do while (abs(psi(N))>tol)
    print*, E, Emin, Emax, psi(N)
    E=(Emin+Emax)/2
    call shoot(x, V, N, E, dx, psi)
    if (psi(N)>0) then
        print*, "first"
        Emin=E
    elseif (psi(N)<0) then 
        print*, "second"
        Emax=E
    endif

enddo

end program bisection    

subroutine potential(x,V,N,omega)
    implicit none
    integer, intent(in) :: N
    real*8, intent(in) :: omega
    real*8, dimension(N), intent(in) :: x
    real*8 , dimension(N), intent(out) :: V
    V=0.5*omega**2*x**2
end subroutine potential

subroutine shoot(x, V, N, E, dx, psi)
    implicit none
    integer :: i, outunit
    real*8, intent(in) :: E, dx
    integer, intent(in)::N
    real*8, dimension(N), intent(in) :: x, V
    real*8, dimension(N), intent(inout) :: psi
    open(newunit=outunit, file='bisection.out', action="write")
    write(outunit,*) x(1), V(1), psi(1), E
    
    do i=2, N, 1
        psi(i+1)=2*(dx**2*(V(i)-E)+1)*psi(i)-psi(i-1)
        write(outunit,*) x(i), V(i), psi(i), E
    enddo
    close(outunit)
end subroutine shoot
