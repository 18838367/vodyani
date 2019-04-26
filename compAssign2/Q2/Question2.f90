program Q2
implicit none
character(len=8) :: conv, f
real*8, allocatable, dimension(:) :: psi1m, psimN, V, x, psi
integer :: N, outunit, i, m, aimNodes, nodes, j
real*8 :: omega, dx, Emax, Emin, xmax, xmin, E, tol, delE
!______CONSTANTS & INITIALISATION
omega=1
dx=0.01
xmax=5
xmin=-5
tol=0.01
nodes=666
N=ceiling((xmax-xmin)/dx)
m=ceiling(N/2.0)-2
allocate(psi1m(m), psimN(N-m), psi(N), V(N), x(N))
do i=1, N
    x(i)=dx*i+xmin
enddo
call potential(x,V,N,omega)
!______CONSTANTS & INITIALISATION


Emax=4
Emin=0
aimNodes=1

!______Checking number of nodes
do while (nodes/=aimNodes)
    E=(Emin+Emax)/2
    call numerovForw(V, psi1m, dx, E, N, m, x, aimNodes)
    call numerovBack(V, psimN, dx, E, N, m, x, aimNodes)
    call countNodes(N, m, psi1m, psimN, nodes)
    print*, E, Emin, Emax
    print*, nodes
    if (nodes>aimNodes) then
        Emax=E
    elseif (nodes<aimNodes) then
        Emin=E
    endif
enddo
!______Checking numer of nodes


!______Cooley Energy Correction

delE=0
do i=1, m
    psi(i)=psi1m(i)
enddo
do i=1, N-m
    psi(m+i)=psimN(i)
enddo
call cooleyE(V, psi, dx, E, N, m, delE)
print*, delE

do while(abs(delE)>tol)
    E=E+abs(delE)
    print*, E
    call numerovForw(V, psi1m, dx, E, N, m, x, aimNodes)
    call numerovBack(V, psimN, dx, E, N, m, x, aimNodes) 
    do i=1, m
        psi(i)=psi1m(i)
    enddo
    do i=1, N-m
        psi(m+i)=psimN(i)
    enddo
    call cooleyE(V, psi, dx, E, N, m, delE)
enddo

!_____Cooley Energy Correction

!_____Format and write to file

f='(I2.2)' !format
write(conv,f) aimNodes !converts N to string

open(newunit=outunit, file='psi'//trim(conv)//'.out', action="write")
do i=1, N
    write(outunit, *) x(i), V(i), psi(i)
enddo
close(outunit)   

!_____Formart and write to file

end program Q2


subroutine numerovForw(V, psi1m, dx, E, N, m, x, aimNodes)
    implicit none
    integer, intent(in) :: N, m, aimNodes
    real*8, dimension(N), intent(in) :: V, x
    real*8, dimension(m), intent(out) :: psi1m
    real*8, intent(in) :: dx, E
    real*8, dimension(N) :: g
    integer :: i, outunit
    psi1m(1)=0
    psi1m(2)=(-1)**(aimNodes)*0.0001

    g=2*(V-E)
    do i=2, m-1
        psi1m(i+1)=(2*(1+5*dx**2/12*g(i))*psi1m(i)-(1-dx**2/12*g(i-1))*psi1m(i-1))/(1-dx**2/12*g(i+1))
    enddo
    open(newunit=outunit, file="numerovForw.out", action="write")
    do i=1, m
        write(outunit,*) x(i), V(i), psi1m(i)
    enddo
    close(outunit)
    print*, size(psi1m), m, N, "------------"
end subroutine numerovForw

subroutine numerovBack(V, psimN, dx, E, N, m, x, aimNodes)
    implicit none
    integer, intent(in) :: N, m, aimNodes
    real*8, dimension(N), intent(in) :: V, x
    real*8, dimension(N-m), intent(out) :: psimN
    real*8, dimension(N-m) :: temp
    real*8, intent(in) :: dx, E
    real*8, dimension(N) :: g
    integer :: i, outunit
    temp(1)=0
    temp(2)=0.0001

    g=2*(V-E)
    do i=2, N-m-1
        temp(i+1)=(2*(1+5*dx**2/12*g(i))*temp(i)-(1-dx**2/12*g(i-1))*temp(i-1))/(1-dx**2/12*g(i+1))
    enddo
    do i=1, N-m
        psimN(N-m-i)=temp(i)
    enddo
    open(newunit=outunit, file="numerovBack.out", action="write")
    do i=1, N-m
        write(outunit,*) x(i+m), V(i+m), psimN(i)
    enddo
    close(outunit)
    print*, size(psimN), m, N, "------------"
end subroutine numerovBack

subroutine potential(x,V,N,omega)
    implicit none
    integer, intent(in) :: N
    real*8, intent(in) :: omega
    real*8, dimension(N), intent(in) :: x
    real*8 , dimension(N), intent(out) :: V
    V=0.5*omega**2*x**2
end subroutine potential


subroutine countNodes(N, m, psi1m, psimN, nodes)
    implicit none
    integer, intent(in) :: N, m
    real*8, dimension(m), intent(in) :: psi1m
    real*8, dimension(N-m), intent(in) :: psimN
    integer, intent(out) :: nodes
    integer :: i
    nodes=0
    do i=2, m
        if(psi1m(i-1)*psi1m(i)<0) then 
            print*, "NODE FOUND:", i, N
            nodes=nodes+1
        endif
    enddo
    do i=2, (N-m)
        if(psimN(i-1)*psimN(i)<0) then
            print*, "NODE FOUND:", i+m, N
            nodes=nodes+1
        endif
    enddo
end subroutine countNodes

subroutine cooleyE(V, psi, dx, E, N, m, delE)
    implicit none
    integer, intent(in) :: N, m
    real*8, intent(in) :: dx, E
    real*8, intent(out) :: delE
    real*8, dimension(N), intent(in) :: psi, V
    real*8, dimension(N) :: Y, g, delEArray

    g=2*(V-E)
    Y=(1-dx**2/12*g)*psi
    delE=psi(m)/sum(psi**2)*(-0.5*(Y(m+1)-2*Y(m)+Y(m-1))/(dx**2)+(V(m)-E)*psi(m))
!    print*, psi(m+50)/sum(psi**2)*(-0.5*(Y(m+50+1)-2*Y(m+50)+Y(m-1+50))/(dx**2)+(V(m+50)-E)*psi(m+50))
!    print*, psi(m+101)/sum(psi**2)*(-0.5*(Y(m+1+101)-2*Y(m+101)+Y(m-1+101))/(dx**2)+(V(m+101)-E)*psi(m+101))
!    print*, psi(m), psi(m+1), psi(m-1)
end subroutine cooleyE
