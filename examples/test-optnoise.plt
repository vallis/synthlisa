#!/usr/local/bin/gnuplot -persist

set terminal x11

set title "Optical-path noise";
set xlabel "f [Hz]"; set ylabel "S(f) [1/Hz, one-sided]";

set logscale x;
set logscale y;

set key left top

set format x "10^{%L}" 
set format y "10^{%L}"

plot "data/optnoise-freq.txt" every 4 using 1:2 title 'Pseudorandom optical-path noise' with lines;
replot 1.8e-37*x**2 title 'Nominal optical-path noise';

set size 0.8,0.8; set ylabel 2,0
set terminal postscript eps enhanced "Times" 18 color;
set output "eps/test-optnoise.eps";
replot;
