program Question3
implicit none 
integer :: N, i, j
real*8, allocatable, dimension(:,:) :: A, B
real*8 :: temp=0, t1, t2
write(*,*) "Please enter the number of dimensions"
read(*,*) N
allocate(A(N,N), B(N,N))
call cpu_time(t1)
do i=1, N, 1
    do j=1, N, 1
        A(i,j)=i*j
    enddo
enddo
call cpu_time(t2)
write(*,*) "Row assignment time =", (t2-t1)

call cpu_time(t1)
do i=1, N, 1
    do j=1, N, 1
        A(j,i)=i*j
    enddo
enddo
call cpu_time(t2)
write(*,*) "Column assignment time =", (t2-t1)

end program Question3
