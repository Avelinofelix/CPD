CC = gcc
CFLAGS = -Wall -Wextra -pedantic -std=c99
LDFLAGS = -lm

TARGET = sisRecom
OBJS = sisRecom.o

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

sisRecom.o: sisRecom.c
	$(CC) $(CFLAGS) -c $^ -o $@

.PHONY: clean run

clean:
	rm -f $(TARGET) $(OBJS)

run: $(TARGET)
	./$(TARGET) $(file) | tee $(file:.in=.out)

