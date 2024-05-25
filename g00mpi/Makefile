CC = mpicc
CFLAGS = -Wall -Wextra -pedantic -std=c99
LDFLAGS = -lm

TARGET = sisRecom-mpi
OBJS = sisRecom-mpi.o

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

sisRecom-mpi.o: sisRecom-mpi.c
	$(CC) $(CFLAGS) -c $^ -o $@

.PHONY: clean run

clean:
	rm -f $(TARGET) $(OBJS)

run: $(TARGET)
	@echo "Digite o número de processos:"
	@read num_procs; \
	mpirun -np $$num_procs ./$(TARGET) $(file) | tee $(file:.in=.out)