#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "Laser frequency noise";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set key left bottom;

set format x "10^{%L}"; 
set format y "10^{%L}";

plot "data/lasernoise-freq.txt" every 4 using 1:2 title 'Pseudorandom laser frequency noise' with lines;

set parametric;

set trange [1e-4:0.5];
replot t,1.1e-26 title 'Nominal laser frequency noise';

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced "Times" 18 color;
set output "eps/test-lasernoise.eps";
replot;
