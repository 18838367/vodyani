set terminal postscript enhanced eps color dashed "Helvetica, 14" size 8cm,7cm
set key bottom right
set output "Q1.eps"
set datafile separator ","
set samp 500
set ylabel "f(x)"
set xlabel "x"

p "p1.out" t "Discrete (Fortran)" with points pt 7 lc 2, log(x)/(1-x) lc 7 t "Continuous (GNUplot)"
