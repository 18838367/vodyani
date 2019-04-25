program bisection
implicit none
integer :: aimNodes
character(len=8) :: conv, f
real*8 :: Emax, Emin, dx
real*8 :: psi1, psi2, E, xmax, xmin, omega, tol
real*8, allocatable, dimension(:) :: x, V, psi
integer :: i, N, outunit, nodes, flag, j

aimNodes=0
dx=0.001
omega=1
flag=0
xmin=-5
xmax=5
N=ceiling((xmax-xmin)/dx)
allocate(x(N), V(N), psi(N))


do j=0, 3
    print*, "jjjjjjj"
    print*, j
    aimNodes=j
    Emax=4
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
    call shoot(x, V, N, E, dx, psi, aimNodes)
    call countNodes(N, psi, nodes)
    print*, E, Emin, Emax, psi(N)
    print*, nodes
    if (nodes>aimNodes) then
        Emax=E
    elseif (nodes<aimNodes) then
        Emin=E
    endif
    
    if(nodes==aimNodes) then
        if (psi(N)>0) then
            print*, "psiHigh", psi(N)
            Emin=E
        elseif (psi(N)<0) then 
            Emax=E
            print*, "psiLow", psi(N)
        endif
    endif
    
    do while (abs(psi(N))>tol)
        E=(Emin+Emax)/2
        call shoot(x, V, N, E, dx, psi, aimNodes)
        call countNodes(N, psi, nodes)
        print*, E, Emin, Emax, psi(N)
        print*, nodes
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
    
    f='(I2.2)' !format
    write(conv,f) aimNodes !converts N to string

    open(newunit=outunit, file='normalised'//trim(conv)//'.out', action="write")
    do i=1, N
        write(outunit,*) x(i), V(i), psi(i)/(sum(psi)*dx)+E
    enddo
    close(outunit)
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

subroutine shoot(x, V, N, E, dx, psi, aimNodes)
    implicit none
    integer :: i, outunit
    real*8, intent(in) :: E, dx
    character(len=8) :: f, conv
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
    do i=2, N
        if(psi(i-1)*psi(i)<0) then 
            nodes=nodes+1
            print*, "NODE FOUND:", i, N, i/N
        endif
    enddo
end subroutine countNodes
