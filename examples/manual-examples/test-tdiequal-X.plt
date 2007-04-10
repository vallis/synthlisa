#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "First-generation TDI noise";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set key left top;

set format x "10^{%L}"; 
set format y "10^{%L}";

plot "data/tdiequal-X-theory.txt"  using 1:2 title 'X, theory' with lines;
replot   "data/tdiequal-X-sampled.txt" using 1:2 title 'X, synthLISA' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced color "Times" 18;
set output "eps/test-tdiequal-X.eps";
replot;
