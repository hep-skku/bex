# Makefile

## Define paths to the libraries
PDFLIB=$(HOME)/work/Blackhole/lhapdf/lhapdf-5.8.9/lib

## Compiler options
CXX=g++
F77=gfortran
EXE=bex
CXXFLAGS=-I.

## Actions
vpath %.h include
vpath %.o tmp

all: $(EXE)

SRCS=$(filter-out main.cc,$(notdir $(wildcard src/*.cc)))
OBJS=$(SRCS:.cc=.o)

$(OBJS): $(addprefix src/,$(*F:.o=.cc))
	$(CXX) $(CXXFLAGS) -o tmp/$@ -c src/$(*F).cc

$(EXE): $(OBJS) src/main.cc
	$(CXX) $(CXXFLAGS) -o bin/$(EXE) src/main.cc $(addprefix tmp/,$(OBJS))

clean:
	-rm -f bin/$(EXE)
	-rm -f tmp/*
