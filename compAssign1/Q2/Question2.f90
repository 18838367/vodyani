program Question2
implicit none 
integer :: N, i, j, k, C=0, numThread, t, maxThreads, outputunit
real*8, allocatable, dimension(:,:) :: A, B, M
real*8 :: temp=0, t1, t2
character(len=8) :: f, conv !format variable and temporary conversion variable
real*8, external :: omp_get_wtime
integer, external :: omp_get_num_threads, omp_get_max_threads

write(*,*) "Please enter the number of dimensions"
read(*,*) N
allocate(A(N,N), B(N,N), M(N,N))

!below code simply converts dimension N into a string for use in the file name
f='(I5.5)' !format
write(conv,f) N !converts N to string

open(newunit=outputunit, file="Q2_N"//trim(conv)//".out", action="write")
write(outputunit,*) "Number of Threads,   Processing Time"
close(outputunit)

!intialising values (1-N^2 across both arrays)
do i=1, N, 1
    do j=1, N, 1
        A(i,j)=1+C
        B(i,j)=1+C
        C=C+1
    enddo
enddo
!write(*,*) A, B

!$OMP PARALLEL
maxThreads=omp_get_max_threads()
!$OMP END PARALLEL
write(*,*) "The maximum number of threads is", maxThreads

do t=1, maxThreads, 1
    call omp_set_num_threads(t)
    !$OMP PARALLEL
    numThread=omp_get_num_threads()
    !$OMP END PARALLEL
    write(*,*) "Processing with ", numThread, " parallel threads"

    !matrix multiplication below
    !################################## TIMED #################

    t1=omp_get_wtime()
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

    !################################## TIMED #################
    !check matrix multiplication
    
    write(*,*) "Matrix multiplcation correct : ", (sum(M/matmul(A,B))==N*N)
    open(newunit=outputunit, file="Q2_N"//trim(conv)//".out", action="write", position="append")
    write(outputunit,*) numThread, (t2-t1)
    close(outputunit)
enddo
end program Question2
