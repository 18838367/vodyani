program Question3
implicit none 
integer :: N, i, j, k, C=0
real*8, allocatable, dimension(:,:) :: A, B
real*8 :: temp=0, t1, t2
write(*,*) "Please enter the number of dimensions"
read(*,*) N
allocate(A(N,N), B(N,N))
call cpu_time(t1)
do i=1, N, 1
    do j=1, N, 1
        A(i,j)=1+C
        C=C+1
    enddo
enddo
call cpu_time(t2)
write(*,*) t2-t1

call cpu_time(t1)
do i=1, N, 1
    do j=1, N, 1
        A(j,i)=1+C
        C=C+1
    enddo
enddo
call cpu_time(t2)
write(*,*) t2-t1

end program Question3
