program Question2
implicit none 
integer :: N, i, j, k, C=0, numThread, t
real*8, allocatable, dimension(:,:) :: A, B, M
real*8 :: temp=0, t1, t2
real*8, external :: omp_get_wtime
integer, external :: omp_get_num_threads

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

do t=1, 4, 1
    call omp_set_num_threads(t)
    !$OMP PARALLEL
    numThread=omp_get_num_threads()
    !$OMP END PARALLEL
    !matrix multiplication below
    t1=omp_get_wtime()
    write(*,*) numThread
    !$OMP PARALLEL DO PRIVATE(temp) 
    do i=1, N, 1
        do j=1, N, 1
            temp=0
            do k=1, N, 1
                temp=temp+(A(i,k)*B(k,j))
            enddo
            M(i,j)=temp
        enddo
    enddo
    !$OMP END PARALLEL DO
    
    t2=omp_get_wtime()
    !check matrix multiplication
    
    write(*,*) (sum(M/matmul(A,B)))
    write(*,*) (t2-t1)
enddo
end program Question2
