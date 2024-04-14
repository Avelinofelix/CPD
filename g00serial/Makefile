CC = gcc
CFLAGS = -Wall -Wextra -pedantic -std=c99
LDFLAGS = -lm

sisRecom: sisRecom.o
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

sisRecom.o: sisRecom.c
	$(CC) $(CFLAGS) -c $<

.PHONY: clean run

clean:
	rm -f sisRecom sisRecom.o

run: sisRecom dado0.in
	./sisRecom dado0.in

