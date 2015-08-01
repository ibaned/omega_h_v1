CFLAGS = -Weverything
sources := $(wildcard src/*.c)
headers := $(wildcard src/*.h)
objects := $(patsubst %.c,%.o,$(sources))
default: $(objects)
clean:
	rm -f *.o
.PHONY: default clean
