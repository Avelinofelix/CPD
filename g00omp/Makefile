CC = gcc
CFLAGS = -Wall -Wextra -pedantic -std=c99 -fopenmp 
LDFLAGS = -lm

TARGET = sisRecom-omp
OBJS = sisRecom-omp.o

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

sisRecom-omp.o: sisRecom-omp.c
	$(CC) $(CFLAGS) -c $^ -o $@

.PHONY: clean run

clean:
	rm -f $(TARGET) $(OBJS)

run: $(TARGET)
	./$(TARGET) $(file) | tee $(file:.in=.out)

