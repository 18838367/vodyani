program Q3
implicit none
character(len=8) :: f, conv
real*8, allocatable :: K(:,:), V(:,:), H(:,:), B(:,:), w(:), z(:, :), r(:), xi(:,:), phi(:,:)
integer :: i, ell, num, N, ierr, j, g
real*8 :: lam, rmax, dr
ell=0
lam=4
read(*,*) num
dr=0.001
rmax=50
N=ceiling(rmax/dr)
allocate(r(N), xi(num,N), phi(num,N))
do i=1, N 
    r(i)=dr*i
enddo
allocate(K(num,num), H(num,num), V(num,num), B(num,num), w(num), z(num, num))
call Kmatrix(num, ell, lam, K)
call Vmatrix(num, ell, lam, V)
call Bmatrix(num, ell, B)
H=K+V
print*, size(H), shape(H)
print*, size(K), shape(K)
print*, size(V), shape(V)
call rsg(num,num,H,B,w,1,z,ierr)
print*, "Eigenvalues :", w
print*, size(z), shape(z)
print*, "Eigenvectors :", z

call LagBasis(lam, ell, num, r, xi, N)
do g=1,num
    do i=1,N
        phi(g,i)=0
    enddo
enddo

do g=1, num
    do i=1,N
        do j=1, num
            phi(g,i)=phi(g,i)+xi(j,i)*z(j,g)
        enddo
    enddo
enddo
f='(I4)' !format
write(conv,f) num !converts to string for use in name

open(unit=500, file="energies"//trim(conv)//".out", action="write")
do j=1, num
    write(500,*) w(j)
enddo
close(500)

open(unit=400, file="phi"//trim(conv)//".out", action="write")
do i=1, N
    write(400,*) r(i), phi(:,i)
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
        if (i==j) then
            B(i,j)=1
        endif
        if ((i-1)==j) then
            B(i,j)=-0.5*((1-(ell*(ell+1)/((i+ell)*(i+ell-1))))**(0.5))
        endif
        if (i==(j-1)) then 
            B(i,j)=-0.5*((1-(ell*(ell+1)/((i+ell)*(i+ell+1))))**(0.5))
        endif
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

integer function del(i,j) result(f)
implicit none
integer, intent(in) :: i, j

if(i==j) then
    f=1
else 
    f=0
endif
end function 

subroutine LagBasis(lam, ell, num, r, xi, N)
implicit none
integer, intent(in) :: ell, num, N
real*8, intent(in) :: r(N), lam
real*8, dimension(num,N), intent(out) :: xi
real*8, dimension(num,N) :: L
real*8, external :: factorial
integer :: i, j
real*8 :: e
e=2.71828

call LagPol(num, 2*ell+1, r*lam, N, L) 
do j=1, num
    do i=1, N 
        xi(j,i)=(((lam*factorial(j-1))/((j+ell)*factorial(j+2*ell)))**(0.5))*(lam*r(i))**(ell+1)*e**(-lam*r(i)*0.5)*L(j,i)
    enddo
enddo

open(unit=100, file="Q1.out", action="write")
do i=1, N
    write(100,*) r(i), xi(:,i)
enddo
close(100)
end subroutine LagBasis

!the inputs, mt1=m-1, a=alpha, x=r, k=size(x), L=return polynomial
subroutine LagPol(n, a, x, k, L)
implicit none
integer, intent(in) :: a, n, k
real*8, dimension(k), intent(in) :: x
real*8, dimension(n, k+1), intent(out) :: L
integer :: i, j

do i=1, k
    L(1, i)=1
enddo


if (n>=1) then 
    do i=1, k
        L(2, i)=1+a-x(i)
    enddo
endif  

if (n>1) then
    do j=2, n-1
        do i=1, k
            L(j+1, i)=(((2.0d0*(j-1)+1.0d0+a-x(i))*L(j, i))-(((j-1)+a)*L(j-1,i)))/((j-1)+1.0d0)
        enddo
    enddo
endif
open(unit=200, file="poly.out", action="write")
do i=1, k
    write(200,*) x(i), L(:,i)
enddo
close(200)


end subroutine LagPol

real*8 function factorial(n) result(f)
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
