lisasim: lisasim-lisa.o lisasim-tdi.o lisasim-tens.o lisasim-wave.o lisasim-main.o lisasim-tdifast.o
	g++ -O3 -Wall lisasim-tdi.o lisasim-lisa.o lisasim-wave.o lisasim-tens.o lisasim-main.o lisasim-tdifast.o -o lisasim

%.o: %.cpp lisasim.h
	g++ -O3 -Wall -c $<

clean:
	rm -rf *.o lisasim

.PHONEY: clean
