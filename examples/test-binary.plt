#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "Monochromatic binary at f = 1.94 mHz, first-generation X";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set key left top;

set format x "10^{%L}"; 
set format y "10^{%L}";

plot   "data/tdibinary-X-circ.txt"  using 1:2 title 'X signal, circular LISA' with lines;
replot "data/tdibinary-X-ecc.txt"   using 1:2 title 'X signal, eccentric LISA' with lines;
replot "data/tdibinary-X-noise.txt" using 1:2 title 'X signal + noise, circular LISA' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced "Times" 18 color;
set output "eps/test-binary.eps";
replot;
