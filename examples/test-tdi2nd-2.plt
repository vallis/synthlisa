set terminal x11

set title "Second-generation TDI laser-noise subtraction for eccentric LISA orbit";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set format x "10^{%L}"; 
set format y "10^{%L}";

plot   "data/tdi2nd-X1.txt"           using 1:2 title 'X1, synthLISA' with lines;
replot "data/tdi2nd-X1-nolaser.txt"   using 1:2 title 'X1, synthLISA, no laser noise' with lines;

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced color;
set output "eps/test-tdi2nd-X1.eps";
replot;
