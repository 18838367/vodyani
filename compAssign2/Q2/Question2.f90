program Question2
implicit none
character(len=8):: conv, f
real*8, allocatable, dimension(:) :: V, x, psiForw, psiBack, psi
real*8 :: xmax, xmin, dx, E, Emax, Emin, omega, delE, tol, norm
integer :: N, i, nodes, aimnodes, m, outunit, j, iterations

xmin=-5
xmax=5
Emin=0
Emax=5
dx=0.00
read(*,*) dx
print*, dx 
omega=1
aimNodes=0
E=0
N=(xmax-xmin)/dx
nodes=666
tol=0.000001

allocate(x(N), V(N), psiForw(N), psiBack(N), psi(N))

do i=1, N
    x(i)=xmin+dx*i
enddo

m=(N/2)+1

call potential(x, V, N, omega)
do j=0, 3
    iterations=0
    aimNodes=j
    Emax=5
    Emin=0
    do while(nodes/=aimNodes)
        E=(Emax+Emin)/2
        call numerovForw(V, psiForw, dx, E, N, x, aimNodes)
        call numerovBack(V, psiBack, dx, E, N, x, aimNodes)
        iterations=iterations+1
        do i=1, N
            if(i<=m) then
                psi(i)=psiForw(i)/psiForw(m)
            else
                psi(i)=psiBack(i)/psiBack(m)
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
    
    !open(newunit=outunit, file="newNoded.out", action="write")
    !do i=1, N
    !    write(outunit,*) x(i), V(i), psi(i), psiForw(i), psiBack(i)
    !enddo
    !close(outunit)
    
    print*, E, delE
    
    call cooleyE(V, psi, dx, E, N, m, delE)
    
    do while(abs(delE)>tol)
        E=E+delE
        print*, E, delE
        call numerovForw(V, psiForw, dx, E, N, x, aimNodes)
        call numerovBack(V, psiBack, dx, E, N, x, aimNodes)
        iterations=iterations+1
        do i=1, N
            if(i<=m) then
                psi(i)=psiForw(i)/psiForw(m)
            else
                psi(i)=psiBack(i)/psiBack(m)
            endif
        enddo
        call cooleyE(V, psi, dx, E, N, m, delE)
    enddo
    f='(8F5.3)' !format
    write(conv,f) aimNodes*dx+dx !converts N to string
    call simpsons(N,dx,psi**2,norm)
    open(unit=j, file="normalised"//trim(conv)//".out", action="write")
    do i=1, N
        if(j==2.or.j==3) then
            write(j,*) x(i), V(i), -psi(i)/(norm**0.5)+E, E, iterations
        else
            write(j,*) x(i), V(i), psi(i)/(norm**0.5)+E, E, iterations
        endif
    enddo
    close(j)
enddo
call analytic()
end program Question2



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

subroutine simpsons(N, dx, f, norm)
    implicit none
    integer, intent(in) :: N
    real*8, intent(in) :: dx, f(N)
    real*8, intent(out) :: norm
    real*8 :: w(N)
    integer :: i
    
    do i=1, N
    w(i) = 2.0d0 + mod(i+1,2)*2.0d0
    enddo
    w(1) = 1.0d0
    w(N) = 1.0d0
    w = w * dx / 3.0d0
    norm=sum(f*w)
end subroutine simpsons

subroutine analytic()
    implicit none
    real*8 :: pi, e, dx, norm
    integer :: i, N
    real*8, allocatable, dimension(:) :: x, psi0, psi1, psi2, psi3
    pi=3.14159
    e=2.71828

    dx=0.001
    N=10/0.001

    allocate(x(N), psi0(N), psi1(N), psi2(N), psi3(N))
    
    do i=1, N
        x(i)=-5+dx*i
    enddo
    
    psi0=pi**(-0.25)*e**(-x**2/2)!+0.5
    psi1=pi**(-0.25)*2**(0.5)*x*e**(-x**2/2)!+1.5
    psi2=pi**(-0.25)*(1/(2**0.5))*(2*x**2-1)*e**(-x**2/2)!+2.5
    psi3=pi**(-0.25)*(1/(3**0.5))*(2*x**3-3*x)*e**(-x**2/2)!+3.5

!    psi0=psi0/(sum(abs(psi0))*dx)
!    psi1=psi1/(sum(abs(psi1))*dx)
!    psi2=psi2/(sum(abs(psi2))*dx)
!    psi3=psi3/(sum(abs(psi3))*dx)

    call simpsons(N,dx,psi0**2, norm)
    print*, norm, "norm"
    call simpsons(N,dx,psi1**2, norm)
    print*, norm, "norm"
    call simpsons(N,dx,psi2**2, norm)
    print*, norm, "norm"
    call simpsons(N,dx,psi3**2, norm)
    print*, norm, "norm"
    psi0=psi0+0.5
    psi1=psi1+1.5
    psi2=psi2+2.5
    psi3=psi3+3.5

    open(unit=787, file="analytic.out", action="write")
    do i=1, N
        write(787,*) x(i), psi0(i), psi1(i), psi2(i), psi3(i)
    enddo
    close(787)
end subroutine analytic
