#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "Laser frequency noise";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set key left bottom;

set format x "10^{%L}"; 
set format y "10^{%L}";

set xrange [1e-4:1e1];
plot 1.1e-26 title 'Nominal laser frequency noise';

replot "data/lasernoise-freq.txt" every 4 using 1:2 title 'Pseudorandom laser frequency noise' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced color "Times" 18;
set output "eps/test-lasernoise.eps";
replot;
