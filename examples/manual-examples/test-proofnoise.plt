#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "Proof-mass noise";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set format x "10^{%L}" 
set format y "10^{%L}"

set xrange [1e-4:1e1];

plot 2.5e-48*x**-2 title 'Nominal proof-mass noise';
replot "data/proofnoise-freq.txt" every 4 using 1:2 title 'Pseudorandom proof-mass noise' with lines;

set size 0.8,0.8; set ylabel 2,0
set terminal postscript eps enhanced color "Times" 18 "Times" 18;
set output "eps/test-proofnoise.eps";
replot;
