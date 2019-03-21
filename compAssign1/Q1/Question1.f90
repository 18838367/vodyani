program Q1
implicit none

integer :: flag=0
real*8 :: x,f

do while(flag==0)
read(*,*) x
    if(x<=0) then
        write(*,*) "INVALID VALUE please try again"
    else
        if(x==1) then 
            write(*,*) "1,-1"
            flag=1
        else
            f=log(x)/(1-x)
            write(*,*) x,",",f
            flag=1
        endif
    endif
enddo
             
end program Q1
