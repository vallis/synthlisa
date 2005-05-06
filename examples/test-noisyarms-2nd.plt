#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "First-generation TDI X noise for noisy armlength measurement (6 at 1e-2)";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set format x "10^{%L}"; 
set format y "10^{%L}";

set key left top;

set xrange [3e-5:1e-1];
set yrange [1e-49:1e-37];

set mxtics 10;

plot   "data/cleanl2-4096.0-30.0--1.txt" using 1:2 title 'X1, 4096.0/30.0/-1.0' with lines;

replot "data/diffl2-4096.0-100.0--1.txt" using 1:2 title 'X1, 4096.0/100.0/-1.0' with lines;
replot "data/diffl2-4096.0-30.0--1.txt" using 1:2 title 'X1, 4096.0/30.0/-1.0' with lines;

replot "data/diffl2-256.0-100.0--1.txt" using 1:2 title 'X1, 256.0/100.0/-1.0' with lines;
replot "data/diffl2-256.0-30.0--1.txt" using 1:2 title 'X1, 256.0/30.0/-1.0' with lines;

#set terminal x11 2

#set logscale x;
#set logscale y;

#set xrange [1e-2:1e-1]
#set yrange [1e0:1e2]

#set key left top;

#plot   "data/bothl2-4096.0-100.0--1.txt" using 1:($2/$3) 
#replot "data/bothl2-4096.0-30.0--1.txt" using 1:($2/$3) 

set size 0.8,0.8; set ylabel 2,0;
set terminal postscript eps enhanced "Times" 18;
#set terminal postscript eps enhanced color "Times" 18;
set output "eps/test-noisyarms-2nd.eps";
replot;
