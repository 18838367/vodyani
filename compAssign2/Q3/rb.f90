program test_RB
implicit none
real*8, external :: riccati_bessel
integer :: ix
real*8 :: x
open(unit=100,file='RB.out',action='write')
do ix=1,4000
x = dble(ix)/100000.0d0
write(100,*) x, riccati_bessel(3,x)
enddo
end program test_RB

real*8 function riccati_bessel(l,x) result(S)
implicit none
integer, intent(in) :: l
real*8, intent(in) :: x
integer, external :: doubleFac
real*8 :: S0, S1
integer :: j

if(x<0.02) then
    S=(x**(l+1))/doubleFac(2*l+1)
else
    S0 = sin(x); S1 = sin(x)/x - cos(x)
    if(l==0) then
    S = S0
    elseif(l==1) then
    S = S1
    else
    do j=2,l
    S = dble(2*j-1)/x * S1 - S0
    S0 = S1; S1 = S
    enddo
    endif
endif
end function riccati_bessel


integer function doubleFac(n) result(f)
implicit none
integer, intent(in) :: n
integer :: i
i=n
f=1
do while(i>1)
    f=f*i
    i=i-2
enddo
end function doubleFac
    
    
