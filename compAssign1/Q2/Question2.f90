program Question2
implicit none 
integer :: N, i, j, k, C=0
real*8, allocatable, dimension(:,:) :: A, B, M
real*8 :: temp=0


write(*,*) "Please enter the number of dimensions"
read(*,*) N
allocate(A(N,N), B(N,N), M(N,N))


!intialising values (1-N^2 across both arrays)

do i=1, N, 1
    do j=1, N, 1
        A(i,j)=1+C
        B(i,j)=1+C
        C=C+1
    enddo
enddo
!write(*,*) A, B
!matrix multiplication below
do i=1, N, 1
    do j=1, N, 1
        temp=0
        do k=1, N, 1
            temp=temp+(A(i,k)*B(k,j))
        enddo
        M(i,j)=temp
    enddo
enddo

!check matrix multiplication

write(*,*) (sum(M/matmul(A,B)))

end program Question2
