program Question5
implicit none
real*8, dimension(3) :: E, B, P, V, F, VxB
integer :: iter, i, outputunit
real*8 :: t, length, step, k
character (LEN=100) :: header

open(unit=100, file="Q5.in", action="read")
read(100,*) header
read(100,*) length, step, E(1), E(2), E(3), B(1), B(2), B(3), P(1), P(2), P(3), V(1), V(2), V(3)

iter=ceiling(length/step)
close(100)

open(newunit=outputunit, file='chargedParticle.out', action="write")
do i=0, iter
    k=0.5*((V(1)**2+V(2)**2+V(3)**2)**(0.5))
    write(outputunit,*) t, P(1), P(2), P(3), V(1), V(2), V(3), k
    t = t + step
    call cross(V,B,VxB)
    F = VxB+E
    P = P + V*step
    V = V + F*step
enddo


end program Question5


subroutine cross(a, b, axb)
    implicit none
    real*8, dimension(3), intent(out) :: axb
    real*8, dimension(3), intent(in) :: a, b

    axb(1) = a(2) * b(3) - a(3) * b(2)
    axb(2) = a(3) * b(1) - a(1) * b(3)
    axb(3) = a(1) * b(2) - a(2) * b(1)
end subroutine cross



