#OPTIONS = -Wall -g 
#OPTIONS = -Wall -O3 -pg -g 
OPTIONS = -Wall -O3 

lisasim: lisasim-lisa.o lisasim-tdi.o lisasim-tens.o lisasim-wave.o lisasim-main.o lisasim-tdifast.o
	g++ $(OPTIONS) lisasim-tdi.o lisasim-lisa.o lisasim-wave.o lisasim-tens.o lisasim-main.o lisasim-tdifast.o -o lisasim

%.o: %.cpp lisasim.h
	g++ $(OPTIONS) -c $<

clean:
	rm -rf *.o lisasim

.PHONEY: clean
