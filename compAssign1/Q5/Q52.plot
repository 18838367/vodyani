set terminal postscript enhanced eps color dashed "Helvetica, 14" size 8cm,7cm
set key bottom right
set output "Q5_01.eps"
set samp 500
set ylabel "y"
set xlabel "x"

p "chargedParticle_01.out" every::1 using 2:3 t "step=0.1" with lines lc 7, "theory.out" every::1 using 2:3 t "Theoretical" with lines lc 3 
