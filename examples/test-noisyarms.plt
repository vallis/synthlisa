#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "First-generation TDI X noise for noisy armlength measurement";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set format x "10^{%L}"; 
set format y "10^{%L}";

set key left top;
#unset key;

set xrange [1e-4:6e-2];

plot   "data/tdinoisy-clean.txt" using 1:2 title 'X, no armlength error' with lines;
replot "data/tdinoisy-noisy.txt" using 1:2 title 'X, Delta L = 50 m (1 sigma)' with lines;
replot "data/tdinoisy-cleaner.txt" using 1:2 title 'X, Delta L = 50 m (3 sigma)' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced color "Times" 18;
set output "eps/test-noisyarms.eps";
replot;
