#!/usr/local/bin/gnuplot -persist

# this plotting script needs the output of both test-noisyarms-2nd.py
# and test-noisyarms-2nd-2.plt

set terminal x11

set title "Second-generation TDI X1 noise for noisy armlength measurement";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set format x "10^{%L}"; 
set format y "10^{%L}";

set xrange [1e-4:6e-2];

set key left top;
#unset key;

plot   "data/tdinoisy-2nd-clean.txt" using 1:2 title 'X1, no armlength error' with lines;
replot "data/tdinoisy-2nd-noisy.txt" using 1:2 title 'X1, Delta L = 50 m (1 sigma)' with lines;
replot "data/tdinoisy-2nd-cleaner.txt" using 1:2 title 'X1, Delta L = 50 m (3 sigma)' with lines;
replot "data/tdinoisy-2nd-cleanerer.txt" using 1:2 title 'X1, Delta L = 15 m (3 sigma)' with lines;
replot "data/tdinoisy-2nd-cleanest.txt" using 1:2 title 'X1, Delta L = 5 m (3 sigma)' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced color "Times" 18;
set output "eps/test-noisyarms-2nd.eps";
replot;
