#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "First-generation TDI laser-noise subtraction for eccentric LISA orbit";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set format x "10^{%L}"; 
set format y "10^{%L}";

set key left top;
#unset key;

set xrange [1e-4:1e-1];
#set yrange [5e-43:5e-38];

plot   "data/tdi2nd-X-nolaser.txt"    using 1:2 title 'X, synthLISA, no laser noise' with lines;
replot "data/tdi2nd-X.txt"            using 1:2 title 'X, synthLISA, nominal laser noise' with lines;
replot "data/tdi2nd-X-lockedless.txt" using 1:2 title 'X, synthLISA, locked, 0.1x rms laser noise' with lines;
replot "data/tdi2nd-X-locked.txt"     using 1:2 title 'X, synthLISA, nominal laser noise, locked' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced "Times" 18;
set output "eps/test-tdi2nd-X-lock.eps";
replot;
