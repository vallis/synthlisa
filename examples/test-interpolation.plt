#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "Spectral noise distortion due to interpolation";
set xlabel "f [Hz]"; set ylabel "S(f) [a.u., one-sided]";

set logscale x;
set logscale y;

set key left bottom;

set format x "10^{%L}"; 
set format y "10^{%L}";

set xrange [5e-3:1e1];
set yrange [1e-10:1e1];

plot   "data/interpnoise0-freq.txt" using 1:2 title 'Nearest-neighbor' with lines;
replot "data/interpnoise1-freq.txt" using 1:2 title 'Linear' with lines;
replot "data/interpnoise2-freq.txt" using 1:2 title 'Lagrange n = 4' with lines;
replot "data/interpnoise4-freq.txt" using 1:2 title 'Lagrange n = 8' with lines;
replot "data/interpnoise8-freq.txt" using 1:2 title 'Lagrange n = 32' with lines;

# a vertical line at the Nyquist frequency needs a gnuplot trick...

set parametric;
set trange [1e-10:99];
replot 0.5,t notitle;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced color "Times" 18 "Times" 18;
set output "eps/test-interpolation.eps";
replot;
