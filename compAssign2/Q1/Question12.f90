program bisection
implicit none
integer :: aimNodes
character(len=8) :: conv, f
real*8 :: Emax, Emin, dx
real*8 :: psi1, psi2, E, xmax, xmin, omega, tol
real*8, allocatable, dimension(:) :: x, V, psi
integer :: i, N, outunit, nodes, flag, j

aimNodes=0

dx=0.000
read(*,*) dx
print*, dx
omega=1
flag=0
xmin=-5
xmax=5
N=ceiling((xmax-xmin)/dx)
print*, N, "NJNNNNNNNNN"
allocate(x(N), V(N), psi(N))


do j=0, 3
    print*, "Node Number"
    print*, j
    aimNodes=j
    Emax=5
    Emin=0
    psi1=0
    psi2=(-1)**aimNodes*0.01
    do i=1, N
        x(i)=dx*(i-1)+xmin
    enddo
    call potential(x,V,N,omega)
    psi(1)=psi1
    psi(2)=psi2
    tol=0.1
    
    
    E=(Emin+Emax)/2
    call shoot(x, V, N, E, dx, psi, aimNodes, j)
    call countNodes(N, psi, nodes)
!    print*, E, Emin, Emax, psi(N)
!    print*, nodes
    if (nodes>aimNodes) then
        Emax=E
    elseif (nodes<aimNodes) then
        Emin=E
    endif
    
    if(nodes==aimNodes) then
        if (psi(N)>0) then
!            print*, "psiHigh", psi(N)
            Emin=E
        elseif (psi(N)<0) then 
            Emax=E
!            print*, "psiLow", psi(N)
        endif
    endif
    
    do while (abs(psi(N))>tol)
        E=(Emin+Emax)/2
        call shoot(x, V, N, E, dx, psi, aimNodes, j)
        call countNodes(N, psi, nodes)
!        print*, E, Emin, Emax, psi(N)
!        print*, nodes
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
    enddo
    
    f='(8F5.3)' !format
    write(conv,f) aimNodes*dx+dx !converts N to string
    open(unit=j, file='normalised'//trim(conv)//'.out', action="write")
    do i=1, N
        write(j,*) x(i), V(i), psi(i)/(sum(abs(psi))*dx)+E
    enddo
    close(j)
enddo
    
call analytic()

end program bisection    

subroutine potential(x,V,N,omega)
    implicit none
    integer, intent(in) :: N
    real*8, intent(in) :: omega
    real*8, dimension(N), intent(in) :: x
    real*8 , dimension(N), intent(out) :: V
    V=0.5*omega**2*x**2
end subroutine potential

subroutine shoot(x, V, N, E, dx, psi, aimNodes, j)
    implicit none
    integer :: i, outunit
    real*8, intent(in) :: E, dx
    character(len=8) :: f, conv
    integer, intent(in)::N, aimNodes, j
    real*8, dimension(N), intent(in) :: x, V
    real*8, dimension(N), intent(inout) :: psi
    !below code simply converts dimension N into a string for use in the file name
    f='(8F5.3)' !format
    write(conv,f) aimNodes*dx+dx !converts N to string
    open(unit=j, file='bisection'//trim(conv)//'.out', action="write")
    write(j,*) x(1), V(1), psi(1), E
    do i=2, N
        psi(i+1)=2*(dx**2*(V(i)-E)+1)*psi(i)-psi(i-1)
        write(j,*) x(i), V(i), psi(i), E
    enddo
    close(j)
end subroutine shoot

subroutine countNodes(N, psi, nodes)
    implicit none
    integer, intent(in) :: N
    real*8, dimension(N), intent(in) :: psi
    integer, intent(out) :: nodes
    integer :: i
    nodes=0
    do i=2, N
        if(psi(i-1)*psi(i)<0) then 
            nodes=nodes+1
!            print*, "NODE FOUND:", i, N, i/N
        endif
    enddo
end subroutine countNodes

subroutine analytic()
    implicit none
    real*8 :: pi, e, dx 
    integer :: outunit, i, N
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

    psi0=psi0/(sum(abs(psi0))*dx)
    psi1=psi1/(sum(abs(psi1))*dx)
    psi2=psi2/(sum(abs(psi2))*dx)
    psi3=psi3/(sum(abs(psi3))*dx)

    psi0=psi0+0.5
    psi1=psi1+1.5
    psi2=psi2+2.5
    psi3=psi3+3.5

    open(newunit=outunit, file="analytic.out", action="write")
    do i=1, N
        write(outunit,*) x(i), psi0(i), psi1(i), psi2(i), psi3(i)
    enddo
    close(outunit)
end subroutine analytic
