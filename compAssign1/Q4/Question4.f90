program fixed_implementation
implicit none
integer :: n, i, Nmax, outunit
real*8 :: exp1, exp2, exp2out, fac, term_n, t1, t2 
real*8, parameter :: x=5.5

open(newunit=outunit, file='good.dat', action='write')
write(outunit,*) 'Nmax, exp(-5.5), 1/exp(5.5), time'
call cpu_time(t1)
exp1=1.0d0 !exp1 = exp(-5.5) 
exp2=1.0d0 !exp2 = 1/exp(5.5) 
fac=1.0d0
call cpu_time(t2)

write(outunit,*) 0, exp1, exp2, (t2-t1)
do Nmax=1,150

    fac=fac*Nmax
    exp1=exp1+((-1)**Nmax)*(x**Nmax)/fac
    exp2=exp2+(x**Nmax)/fac
    exp2out = 1.0d0/exp2 

    call cpu_time(t2)
    write(outunit,*) Nmax, exp1, exp2out, (t2-t1) 
enddo
close(outunit)

end program fixed_implementation
