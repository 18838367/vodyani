set terminal postscript enhanced eps color dashed "Helvetica, 14" size 8cm,7cm
set key top right
set output "3D.eps"
set samp 500
set ylabel "y"
set xlabel "x"
set zlabel "z"
set view 70,80
splot "chargedParticle.out" every ::1 using 2:3:4 with lines t "Charged particle trajectory"
