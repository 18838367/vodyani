program Q1
implicit none
real*8, allocatable :: r(:), xi(:,:)
integer :: i, ell, num, N
real*8 :: lam, rmax, dr
dr=0.1
ell=1
rmax=10+ell
N=ceiling(rmax/dr)
num=6
lam=2
allocate(r(N), xi(num,N))
do i=1, N 
    r(i)=dr*i
enddo
call LagBasis(lam,ell,num,r,xi,N) 
end program


subroutine LagBasis(lam, ell, num, r, xi, N)
implicit none
integer, intent(in) :: ell, num, N
real*8, intent(in) :: r(N), lam
real*8, dimension(num,N), intent(out) :: xi
real*8, dimension(num,N) :: L
integer, external :: factorial
integer :: i, j
real*8 :: e

e=2.71828

call LagPol(num-1, 2*ell+1, r, N, L) 
do j=1, num
    do i=1, N
        xi(j,i)=((lam*factorial(j-1))/((j+1)*factorial(j+2*ell)))**(0.5)*(lam*r(i))**(ell+1)*e**(-lam*r(i)*0.5)*L(j,i)
    enddo
enddo

open(unit=100, file="Q1.out", action="write")
do i=1, N
    write(100,*) xi(:,i)
enddo
close(100)
end subroutine LagBasis

!the inputs, mt1=m-1, a=alpha, x=r, k=size(x), L=return polynomial
subroutine LagPol(mt1, a, x, k, L)
implicit none
integer, intent(in) :: a, mt1, k
real*8, dimension(k), intent(in) :: x
real*8, dimension(mt1+1, k), intent(out) :: L
integer :: i, j


do i=1, k
    L(1, i)=1
enddo

if (mt1>=1) then 
    do i=1, k
        L(2, i)=1+a-x(i)
    enddo
endif  

if (mt1>1) then
    do j=2, mt1-1
        do i=1, k
            L(j+1, i)=((2*j+1+a-x(i))*(L(j, i))-(j+a)*L(j-1,i))/(j+1)
        enddo
    enddo
endif
open(unit=200, file="poly.out", action="write")
do i=1, k
    write(200,*) L(:,i)
enddo
close(100)


end subroutine LagPol

integer function factorial(n) result(f)
implicit none
integer, intent(in) :: n
integer :: i
i=n
f=1
do while(i>1)
    f=f*i
    i=i-1
enddo
end function factorial


