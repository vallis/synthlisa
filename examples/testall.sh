#!/bin/sh

for file in test-*.py; do python $file; done
for file in test-*.plt; do gnuplot -persist $file; done
