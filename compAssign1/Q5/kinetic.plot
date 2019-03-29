set terminal postscript enhanced eps color dashed "Helvetica, 14" size 8cm,7cm
set key top right
set output "kinetic0.1.eps"
set samp 500
set ylabel "Kinetic Energy"
set xlabel "t"

p "chargedParticle_01.out" every ::1 using 1:8 with lines lc 1
set output "kinetic0.01.eps"
p "chargedParticle_001.out" every ::1 using 1:8 with lines lc 1
set output "kinetic0.00001.eps"
set yrange [0:1]
p "chargedParticle_00001.out" every ::1 using 1:8 with lines lc 1

