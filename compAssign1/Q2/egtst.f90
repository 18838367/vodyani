program OMP_example_1
implicit none 
real*8, allocatable :: A(:,:)
integer :: n, m
real*8 :: t1, t2
real*8, external :: omp_get_wtime

allocate(A(40000,40000)) 
t1=omp_get_wtime()
!$OMP PARALLEL DO
    do m=1, 40000 
        do n=m, 40000
            A(n,m) = exp(dble(n*m))
            A(m,n) = A(n,m)
        enddo
 enddo
!$OMP END PARALLEL DO
t2=omp_get_wtime()

print*, 'time taken:', t2-t1
end program
