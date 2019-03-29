program fixed_implementation
implicit none
integer :: n, i, Nmax, outunit
real*8 :: exp1, exp2, fac, term_n, t1, t2 
real*8, parameter :: x=5.5
open(newunit=outunit, file='bad.dat', action='write')
write(outunit,*) 'Nmax, exp(-5.5), 1/exp(5.5), time'
call cpu_time(t1)

do Nmax=0,150
    exp1=1.0d0 !exp1 = exp(-5.5) 
    exp2=1.0d0 !exp2 = 1/exp(5.5) 
    do n=1,Nmax
        fac=1.0d0 
        do i=1,n
            fac=fac*dble(i) 
        enddo

        if(mod(n,2)==0) then 
            exp1 = exp1 + x**n/fac
        else
            exp1 = exp1 - x**n/fac
        endif

        exp2 = exp2 + x**n/fac

    enddo
    exp2 = 1.0d0/exp2 

    call cpu_time(t2)
    write(outunit,*) Nmax, exp1, exp2, (t2-t1) 
enddo
close(outunit)

end program fixed_implementation
