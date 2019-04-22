program Q1
implicit none 
integer :: N0, N1, N2, N3
real*8 :: dx, Emax, Emin
real*8, allocatable, dimension(:) :: psi0, psi1, psi2, psi3

N0=1
N1=3
N2=5
N3=7

dx=0.1
Emin=0
Emax=4

call(bisection(dx, Emax, Emin, N0, psi0)
call(bisection(dx, Emax, Emin, N1, psi1)
call(bisection(dx, Emax, Emin, N2, psi2)
call(bisection(dx, Emax, Emin, N3, psi3)

end program Q1

subroutine bisection(dx, Emax, Emin, aimNodes, psi)
implicit none
real*8, intent(in) :: Emax, Emin, dx
real*8 allocatable, dimension(:), intent(out) :: psi
real*8 :: psi1, psi2, E, xmax, xmin, omega, tol
character(len=8) :: f, conv
real*8, allocatable, dimension(:) :: x, psi, V
integer :: i, N, outunit, nodes, flag
omega=1
flag=0
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
call shoot(x, V, N, E, dx, psi, aimNodes)
call countNodes(N, psi, nodes)
if (nodes>aimNodes) then
    Emax=E
elseif (nodes<aimNodes) then
    Emin=E
endif

if(nodes==aimNodes) then
    if (psi(N)>0) then
        Emin=E
    elseif (psi(N)<0) then 
        Emax=E
    endif
endif

do while (abs(psi(N))>tol)
    print*, E, Emin, Emax, psi(N)
    E=(Emin+Emax)/2
    call shoot(x, V, N, E, dx, psi, aimNodes)
    call countNodes(N, psi, nodes)
    if (psi(N)>0) then
        print*, "first"
        Emin=E
    elseif (psi(N)<0) then 
        print*, "second"
        Emax=E
    endif

enddo

end subroutine bisection    

subroutine potential(x,V,N,omega)
    implicit none
    integer, intent(in) :: N
    real*8, intent(in) :: omega
    real*8, dimension(N), intent(in) :: x
    real*8 , dimension(N), intent(out) :: V
    V=0.5*omega**2*x**2
end subroutine potential

subroutine shoot(x, V, N, E, dx, psi, aimNodes)
    implicit none
    integer :: i, outunit
    real*8, intent(in) :: E, dx
    integer, intent(in)::N, aimNodes
    real*8, dimension(N), intent(in) :: x, V
    real*8, dimension(N), intent(inout) :: psi
    !below code simply converts dimension N into a string for use in the file name
    f='(I2.2)' !format
    write(conv,f) aimNodes !converts N to string

    open(newunit=outunit, file='bisection'//trim(conv)//'.out', action="write")
    write(outunit,*) x(1), V(1), psi(1), E
    
    do i=2, N, 1
        psi(i+1)=2*(dx**2*(V(i)-E)+1)*psi(i)-psi(i-1)
        write(outunit,*) x(i), V(i), psi(i), E
    enddo
    close(outunit)
end subroutine shoot

subroutine countNodes(N, psi, nodes)
    implicit none
    integer, intent(in) :: N
    real*8, dimension(N), intent(in) :: psi
    integer, intent(out) :: nodes
    integer :: i
    nodes=0
    do i=3, N
        if(abs(psi(i-1))>abs(psi(i))) then 
            if(abs(psi(i-1))>abs(psi(i-2))) then
                nodes=nodes+1
            endif
        endif
    enddo
end subroutine countNodes
