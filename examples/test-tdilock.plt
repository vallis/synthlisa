#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "First-generation TDI noise";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set key left top;

plot   "tdi-X-nonlocked.txt"  using 1:2 title 'X, nonlocked' with lines;
replot "tdi-X-locked.txt" using 1:2 title 'X, locked' with lines;
replot "tdi-X-locksimp.txt"  using 1:2 title 'X, simplified, locked' with lines;
# replot "tdi-y-zero.txt" using 1:2 titl 'y_231, for lock 2' with lines;

set terminal postscript color;
set output "test-tdilock-X.eps";
replot;
