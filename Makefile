# Makefile

## Define paths to the libraries
PDFLIB=$(HOME)/work/Blackhole/lhapdf/lhapdf-5.8.9/lib
PYTHIA8=$(HOME)/work/Blackhole/pythia8

## Compiler options
CC=g++
F77=gfortran
EXE=bex
CCFLAGS=-I.
LDFLAGS=-L$(PDFLIB) -lm

## Actions
all: $(EXE)

SRCS=$(filter-out main.cc,$(notdir $(wildcard src/*.cc)))
OBJS=$(SRCS:.cc=.o)

$(OBJS): $(addprefix src/,$(*F:.o=.cc))
	$(CC) $(CCFLAGS) -o tmp/$@ -c src/$(*F).cc

$(EXE): $(OBJS) src/main.cc
	$(CC) $(CCFLAGS) $(LDFLAGS) -o bin/$(EXE) src/main.cc $(addprefix tmp/,$(OBJS))

clean:
	-rm -f bin/$(EXE)
	-rm -f tmp/*
