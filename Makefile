# Makefile

## Define paths to the libraries
LHAPDF=$(HOME)/work/Blackhole/lhapdf/lhapdf-5.8.9
PYTHIA8=$(HOME)/work/Blackhole/pythia8
#INCLUDES=. /opt/local/include $(LHAPDF)/include
INCLUDES=. /opt/local/include

## Compiler options
#CC=g++
CC=clang++
EXE=bex
CCFLAGS=$(addprefix -I,$(INCLUDES))
LDFLAGS=-L$(LHAPDF)/lib -lLHAPDF -lm

## Actions
all: $(EXE)

SRCS=$(filter-out main.cc PDFInterface.cc,$(notdir $(wildcard src/*.cc)))
OBJS=$(addprefix tmp/,$(SRCS:.cc=.o))

$(OBJS):
	$(CC) $(CCFLAGS) -o $@ -c src/$(*F).cc

tmp/PDFInterface.o:
	$(CC) $(CCFLAGS) -I$(LHAPDF)/include -o $@ -c src/$(*F).cc

$(EXE): $(OBJS) tmp/PDFInterface.o
	$(CC) $(CCFLAGS) $(LDFLAGS) -o bin/$(EXE) src/main.cc $(OBJS)

clean:
	-rm -f bin/$(EXE)
	-rm -f tmp/*

