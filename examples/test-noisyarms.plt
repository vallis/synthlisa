#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "First-generation TDI X noise for noisy armlength measurement";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set key left top;

plot   "data/tdinoisy-clean.txt" using 1:2 title 'X, clean' with lines;
replot "data/tdinoisy-noisy.txt" using 1:2 title 'X, noisy' with lines;
replot "data/tdinoisy-cleaner.txt" using 1:2 title 'X, cleaner' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced color;
set output "eps/test-noisyarms.eps";
replot;
