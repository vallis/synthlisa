#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "First-generation TDI noise, arms = (15,16,17) s";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set key left top;

set format x "10^{%L}"; 
set format y "10^{%L}";

set xrange [1e-5:2e-1];

plot   "data/tdiunequal-X-freq.txt"     using 1:2 title 'X'     with lines;
replot "data/tdiunequal-alpha-freq.txt" using 1:2 title 'alpha' with lines;
replot "data/tdiunequal-zeta-freq.txt"  using 1:2 title 'zeta'  with lines;
replot "data/tdiunequal-U-freq.txt"     using 1:2 title 'U'  with lines;

# we're not showing this since it's way too high...
# replot "data/tdiunequal-y231-freq.txt"  using 1:2 title 'y_321' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced "Times" 18 color;
set output "eps/test-tdiunequal.eps";
replot;
