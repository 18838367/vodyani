program tst
implicit none
integer :: n, f
integer :: i
read(*,*) n

i=n
f=1
do while(i>1)
    f=f*i
    i=i-2
enddo
print*, f
end program tst
