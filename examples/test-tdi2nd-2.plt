set terminal x11

set title "Second-generation TDI laser-noise subtraction for eccentric LISA orbit";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set format x "10^{%L}"; 
set format y "10^{%L}";

#unset key;

set xrange [1e-4:1e-1];
set yrange [1e-54:1e-36];

plot   "data/tdi2nd-X1-nolaser.txt"  using 1:2 title 'X1, synthLISA, no laser noise' with lines;
#replot "data/tdi2nd-X1.txt"          using 1:2 title 'X1, synthLISA, nominal laser noise' with lines;
replot "data/tdi2nd-X1-residual.txt" using 1:2 title 'X1, residual laser noise' with lines;

set size 0.8,0.8; set ylabel 2,0;
#set terminal postscript eps enhanced color "Times" 18;
set terminal postscript eps enhanced "Times" 18;
set output "eps/test-tdi2nd-z1.eps";
replot;
