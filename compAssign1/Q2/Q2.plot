set terminal postscript enhanced eps color dashed "Helvetica, 14" size 8cm,7cm
set key top right
set output "Q2_N00100.eps"
set ylabel "Processing Time"
set xlabel "Number of Threads"

p "Q2_N00100.out" every ::1 using 1:2 with linespoints pt 7 lc 7 
set output "Q2_N01000.eps"
p "Q2_N01000.out" every ::1 using 1:2 with linespoints pt 7 lc 7

