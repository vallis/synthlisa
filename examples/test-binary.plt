#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "Monochromatic binary at f = 1.944 mHz, first-generation X";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set key left top;

set format x "10^{%L}"; 
set format y "10^{%L}";

set xrange [0.00194275:0.0019455];
set yrange [1e-43:1e-39];

plot   "data/tdibinary-X-circ.txt"  using 1:2 title 'X, 2 * signal spectral density' with lines;
replot "data/tdibinary-X-noise.txt" using 1:2 title 'X, noise spectral density' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced color "Times" 18 "Times" 18;
set output "eps/test-binary.eps";
replot;
