set terminal postscript enhanced eps color dashed "Helvetica, 14" size 8cm,7cm
set key bottom right
set output "Q0001.eps"
set samp 500
set ylabel "f(x)"
set xlabel "x"


plot "normalised0.001.out" u 1:3 with points lc 1, "normalised0.002.out" u 1:3 with points lc 2, "normalised0.003.out" u 1:3 with points lc 3, "normalised0.004.out" u 1:3 with points lc 4, "analytic.out" using 1:2 with lines lc 8, "analytic.out" using 1:3 with lines lc 8, "analytic.out" using 1:4 with lines lc 8, "analytic.out" using 1:5 with lines lc 8
set output "Q0010.eps"
plot "normalised0.010.out" u 1:3 with points lc 1, "normalised0.020.out" u 1:3 with points lc 2, "normalised0.030.out" u 1:3 with points lc 3, "normalised0.040.out" u 1:3 with points lc 4, "analytic.out" using 1:2 with lines lc 8, "analytic.out" using 1:3 with lines lc 8, "analytic.out" using 1:4 with lines lc 8, "analytic.out" using 1:5 with lines lc 8
set output "Q0100.eps"
plot "normalised0.100.out" u 1:3 with points lc 1, "normalised0.200.out" u 1:3 with points lc 2, "normalised0.300.out" u 1:3 with points lc 3, "normalised0.400.out" u 1:3 with points lc 4, "analytic.out" using 1:2 with lines lc 8, "analytic.out" using 1:3 with lines lc 8, "analytic.out" using 1:4 with lines lc 8, "analytic.out" using 1:5 with lines lc 8
set output "Q1000.eps"
plot "normalised1.000.out" u 1:3 with points lc 1, "normalised2.000.out" u 1:3 with points lc 2, "normalised3.000.out" u 1:3 with points lc 3, "normalised4.000.out" u 1:3 with points lc 4, "analytic.out" using 1:2 with lines lc 8, "analytic.out" using 1:3 with lines lc 8, "analytic.out" using 1:4 with lines lc 8, "analytic.out" using 1:5 with lines lc 8
