###############################################################################
# Makefile for PROJECT
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -O3 -funroll-loops 
LIBS = -lm
FIG = histogram.png
RESU = results.txt

BIN = main

all: $(BIN)

main: main.c
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)
	
clean:
	$(RM) $(BIN) $(FIG) $(RESU)
