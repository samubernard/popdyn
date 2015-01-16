CFLAGS=-Wall -g `pkg-config --cflags gsl` 
LDFLAGS=-g `pkg-config --libs gsl` -lm


all: coupled_oscillators.out

coupled_oscillators.out: coupled_oscillators.c
	gcc $(LDFLAGS) $(CFLAGS) -O -o coupled_oscillators.out  coupled_oscillators.c

clean:
	rm -f *.o

veryclean:
	rm -f *.o; rm -f *.out; rm -rf *.out.dSYM
