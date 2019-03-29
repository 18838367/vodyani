set terminal postscript enhanced eps color dashed "Helvetica, 14" size 8cm,7cm
set key top right
set samp 500
set ylabel "f(x)"
set xlabel "N"
set output "Q4T.eps"
p "bad.dat" every ::1 using 4 t "Time complexity of bad implementation" with points pt 7 lc 2, "good.dat" every::1 using 4 t "Time complexity of my implementation" with points pt 7 lc 3
set xrange [1:20]
set output "Q41.eps"
p "good.dat" every ::1 using 2 t "approximation (Fortran)" with points pt 7 lc 3, 1/(2.71828**(5.5)) lc 7 with lines t "Analytic value (GNUplot)"
set output "Q42.eps"
p "good.dat" every ::1 using 3 t "approximation (Fortran)" with points pt 7 lc 3, 1/(2.71828**(5.5)) lc 7 with lines t "Analytic value (GNUplot)"
