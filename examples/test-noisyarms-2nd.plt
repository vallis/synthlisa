#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "Second-generation TDI X1 noise for noisy armlength measurement, eccentric LISA orbits";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set key left top;

plot   "data/tdinoisy-2nd-clean.txt" using 1:2 title 'X1, clean' with lines;
replot "data/tdinoisy-2nd-noisy.txt" using 1:2 title 'X1, noisy' with lines;
replot "data/tdinoisy-2nd-cleaner.txt" using 1:2 title 'X1, cleaner' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced color;
set output "eps/test-noisyarms-2nd.eps";
replot;
