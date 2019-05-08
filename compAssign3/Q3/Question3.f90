program Q3
implicit none
real*8, allocatable :: K(:,:), V(:,:), H(:,:), B(:,:), w(:), z(:), r(:), xi(:,:), phi(:)
integer :: i, ell, num, N, ierr, j
real*8 :: lam, rmax, dr
ell=1
lam=2
num=5
dr=0.1
rmax=10
N=ceiling(rmax/dr)
allocate(r(N), xi(num,N), phi(N))
do i=1, N 
    r(i)=dr*i
enddo
allocate(K(num,num), H(num,num), V(num,num), B(num,num), w(num), z(num))
call Kmatrix(num, ell, lam, K)
call Vmatrix(num, ell, lam, V)
call Bmatrix(num, ell, B)
H=K+V

call rsg(num,num,H,B,w,1,z,ierr)
print*, "Eigenvalues :", w
print*, "Eigenvectors :", z

call LagBasis(lam, ell, num, r, xi, N)
do i=1,N
    phi(i)=0
enddo
open(unit=400, file="phi.out", action="write")
do i=1,N
    do j=1, num
        phi(i)=phi(i)+xi(j,i)
    enddo
    write(400,*) phi(i) 
enddo
close(400)
print*, "Done!"

end program Q3


subroutine Bmatrix(nm,ell, B)
implicit none
integer, intent(in) :: ell, nm
real*8, intent(out) :: B(nm, nm)
integer :: i, j
integer, external :: del
real*8 :: temp1, temp2
do i=1, nm
    do j=1, nm
        temp1=del(i,j)-0.5*((1-(ell*(ell+1)/((i+ell)*(i+ell+1))))**(0.5))*del(i,j-1)
        temp2=-0.5*((1-(ell*(ell+1)/((i+ell)*(i+ell-1))))**(0.5))*del(i-1,j)
        B(i,j)=temp1+temp2
    enddo
enddo

end subroutine Bmatrix

subroutine Kmatrix(nm, ell, lam, K)
implicit none
integer, intent(in) :: ell, nm
real*8, intent(in) :: lam
real*8, intent(out) :: K(nm,nm)
integer :: i, j
real*8 :: B(nm,nm)
integer, external :: del

call Bmatrix(nm, ell, B)

do i=1, nm
    do j=1, nm
        K(i,j)=-lam**2/8*B(i,j)+lam**2/4*del(i,j)
    enddo
enddo

end subroutine Kmatrix

subroutine Vmatrix(nm, ell, lam, V)
implicit none
integer, intent(in) :: ell, nm
real*8, intent(in) :: lam
real*8, intent(out) :: V(nm, nm)
integer :: i, j
integer, external :: del

do i=1, nm
    do j=1, nm
        V(i,j)=-lam/(2*(i+ell))*del(i,j)
    enddo
enddo

end subroutine Vmatrix

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

integer function del(i,j) result(f)
implicit none
integer, intent(in) :: i, j

if(i==j) then
    f=1
else 
    f=0
endif
end function 
