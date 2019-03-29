set terminal postscript enhanced eps color dashed "Helvetica, 14" size 8cm,7cm
set key bottom right
set output "Q5.eps"
set samp 500
set ylabel "f(x)"
set xlabel "x"

p "chargedParticle_01.out" every::1 using 2:3 t "step=0.1" with lines lc 7, "chargedParticle_001.out" every::1 using 2:3 t "step=0.01" with lines lc 1, "chargedParticle_00001.out" every::1 using 2:3 t "step=0.0001" with lines lc 9 ,"theory.out" every::1 using 2:3 t "Theoretical" with lines lc 3 
