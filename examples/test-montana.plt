#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "Monochromatic binary at f = 1.944 mHz, first-generation X";
set xlabel "f [mHz]"; set ylabel "S(f) [1/Hz, one-sided]";

# set logscale x;
set logscale y;

set key left top;

set format x "%.4t"; 
set format y "10^{%L}";

set xrange [0.00194275:0.0019455];

plot   'data/tdibinary-X-circ-closeup.txt'    using 1:2 title 'X, Synthetic LISA' with lines;
replot 'data/tdibinary-X-montana-closeup.txt' using 1:2 title 'X, LISA Simulator v. 2.0' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced color "Times" 18 "Times" 18;
set output "eps/test-montana.eps";
replot;
