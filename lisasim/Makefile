lisasim: lisasim-lisa.o lisasim-tdi.o lisasim-tens.o lisasim-wave.o lisasim-main.o
	g++ -O3 lisasim-tdi.o lisasim-lisa.o lisasim-wave.o lisasim-tens.o lisasim-main.o -o lisasim

%.o: %.cpp lisasim.h
	g++ -O3 -c $<

