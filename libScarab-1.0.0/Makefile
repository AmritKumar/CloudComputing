CC       = gcc
CFLAGS   = -Wall -g3 -ggdb -std=c99 -I/opt/local/include
LDFLAGS += -lgmp -lflint
SOURCE   = $(shell find * -name '*.c')
OBJECTS  = $(SOURCE:.c=.o)
.PHONY:    clean

all: integer-fhe

integer-fhe: $(OBJECTS)
	$(CC) -o integer-fhe $^ $(LDFLAGS)

clean:
	rm -f $(OBJECTS) integer-fhe

