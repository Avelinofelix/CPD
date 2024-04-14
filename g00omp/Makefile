CC = gcc
CFLAGS = -Wall -Wextra -pedantic -std=c99 -fopenmp
LDFLAGS = -lm -fopenmp

sisRecom-omp: sisRecom-omp.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

sisRecom-omp.o: sisRecom-omp.c
	$(CC) $(CFLAGS) -c $<

.PHONY: clean run

clean:
	rm -f sisRecom-omp sisRecom-omp.o

run: sisRecom-omp dado0.in
	./sisRecom-omp dado0.in
