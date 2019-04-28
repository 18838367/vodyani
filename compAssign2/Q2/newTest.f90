program new
implicit none
real*8, allocatable, dimension(:) :: V, x, psiForw, psiBack, psi
real*8 :: xmax, xmin, dx, E, Emax, Emin, omega, delE, tol
integer :: N, i, nodes, aimnodes, m, outunit

xmin=-5
xmax=5
Emin=0
Emax=5
dx=0.001
omega=1
aimNodes=3
E=0
N=(xmax-xmin)/dx
nodes=666
tol=0.001

allocate(x(N), V(N), psiForw(N), psiBack(N), psi(N))

do i=1, N
    x(i)=xmin+dx*i
enddo

m=(N/2)+1

call potential(x, V, N, omega)

do while(nodes/=aimNodes)
    E=(Emax+Emin)/2
    call numerovForw(V, psiForw, dx, E, N, x, aimNodes)
    call numerovBack(V, psiBack, dx, E, N, x, aimNodes)
    do i=1, N
        if(i<=m) then
            psi(i)=psiForw(i)
        else
            psi(i)=psiBack(i)
        endif
    enddo
    call countNodesAlt(N, m, psi, nodes)
    print*, nodes
    print*, E
    if(nodes>aimNodes) then
        Emax=E
        print*, "top"
    endif
    if(nodes<aimNodes) then
        print*, "bottom"
        Emin=E
    endif
enddo

open(newunit=outunit, file="newNoded.out", action="write")
do i=1, N
    write(outunit,*) x(i), V(i), psi(i), psiForw(i), psiBack(i)
enddo
close(outunit)

print*, E, delE

call cooleyE(V, psi, dx, E, N, m, delE)

do while(abs(delE)>tol)
    E=E+delE
    print*, E, delE
    call numerovForw(V, psiForw, dx, E, N, x, aimNodes)
    call numerovBack(V, psiBack, dx, E, N, x, aimNodes)
    do i=1, N
        if(i<=m) then
            psi(i)=psiForw(i)
        else
            psi(i)=psiBack(i)
        endif
    enddo
    call cooleyE(V, psi, dx, E, N, m, delE)
enddo

open(newunit=outunit, file="newPsi.out", action="write")
do i=1, N
    write(outunit,*) x(i), V(i), psi(i)
enddo
close(outunit)

end program



subroutine numerovForw(V, psi, dx, E, N, x, aimNodes)
    implicit none
    integer, intent(in) :: N, aimNodes
    real*8, dimension(N), intent(in) :: V, x
    real*8, dimension(N), intent(out) :: psi
    real*8, intent(in) :: dx, E
    real*8, dimension(N) :: g
    integer :: i
    psi(1)=0
    psi(2)=(-1)**(aimNodes)*0.0001
    g=2*(V-E)
    do i=2, N-1
        psi(i+1)=(2*(1+5*dx**2/12*g(i))*psi(i)-(1-dx**2/12*g(i-1))*psi(i-1))/(1-dx**2/12*g(i+1))
    enddo
end subroutine numerovForw

subroutine numerovBack(V, psi, dx, E, N, x, aimNodes)
    implicit none
    integer, intent(in) :: N, aimNodes
    real*8, dimension(N), intent(in) :: V, x
    real*8, dimension(N), intent(inout) :: psi
    real*8, intent(in) :: dx, E
    real*8, dimension(N) :: g
    integer :: i
    psi(N)=0
    psi(N-1)=0.0001
    g=2*(V-E)
    do i=1, N-2
        psi(N-i-1)=(psi(N-i+1)*(1-dx**2/12*g(N-i+1))-2*(1+5*dx**2/12*g(N-i))*psi(N-i))/(-1*(1-dx**2/12*g(N-i-1)))
       ! psimN(N-m-i)=(psimN(N-m-i+1)*(1-dx**2/12*g(N-i+1))-2*(1+5*dx**2/12)*psimN(N-m-i))/(-1*(1-dx**2/12*g(N-i-1)))
    enddo
end subroutine numerovBack

subroutine potential(x,V,N,omega)
    implicit none
    integer, intent(in) :: N
    real*8, intent(in) :: omega
    real*8, dimension(N), intent(in) :: x
    real*8 , dimension(N), intent(out) :: V
    V=0.5*omega**2*x**2
end subroutine potential

subroutine countNodes(N, psi, nodes)
    implicit none
    integer, intent(in) :: N
    real*8, dimension(N), intent(in) :: psi
    integer, intent(out) :: nodes
    integer :: i, outunit
    nodes=0
    do i=2, N
        if((psi(i-1)*psi(i))<0) then 
            nodes=nodes+1
        endif
    enddo
end subroutine countNodes

subroutine countNodesAlt(N, m, psi, nodes)
    implicit none
    integer, intent(in) :: N, m
    real*8, dimension(N), intent(in) :: psi
    integer, intent(out) :: nodes
    integer :: i, outunit
    nodes=0
    do i=2, m
        if((psi(i-1)*psi(i))<0) then 
            nodes=nodes+1
        endif
    enddo
    do i=m+2, N
        if((psi(i-1)*psi(i))<0) then
            nodes=nodes+1
        endif
    enddo
end subroutine countNodesAlt

subroutine cooleyE(V, psi, dx, E, N, m, delE)
    implicit none
    integer, intent(in) :: N, m
    real*8, intent(in) :: dx, E
    real*8, intent(out) :: delE
    real*8, dimension(N), intent(in) :: psi, V
    real*8, dimension(N) :: Y, g
    integer :: outunit, i
!    print*, E, "the big one"
    g=2*(V-E)
    Y=(1-dx**2/12*g)*psi
    delE=psi(m)/sum((psi)**2)*(-0.5*(Y(m+1)-2*Y(m)+Y(m-1))/(dx**2)+(V(m)-E)*psi(m))
!    print*, "_____ within"
!    print*, delE, psi(m), psi(m-1), psi(m+1), sum(psi**2)
!    print*, V(m), V(m-1), V(m+1)
!    print*, Y(m), Y(m-1), Y(m+1)
!    print*, g(m), g(m-1), g(m+1)
end subroutine cooleyE
