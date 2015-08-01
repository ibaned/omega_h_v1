CFLAGS = -Weverything
sources := $(wildcard *.c)
objects := $(patsubst %.c,%.o,$(sources))
default: $(objects)
clean:
	rm -f *.o
.PHONY: default clean
