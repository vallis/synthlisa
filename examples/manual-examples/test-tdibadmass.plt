#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "First-generation TDI noise degradation with bad proof mass";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set key left top;

set format x "10^{%L}"; 
set format y "10^{%L}";

set xrange [5e-4:1.25e-1]

plot   "data/tdibadmass-bad.txt" using 1:2 title 'X, bad proof mass 1' with lines;
replot "data/tdibadmass-bad2.txt"  using 1:2 title 'Y, bad proof mass 1' with lines;
replot "data/tdibadmass-bad3.txt"  using 1:2 title 'Z, bad proof mass 1' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced color "Times" 18;
set output "eps/test-tdibadmass.eps";
replot;
