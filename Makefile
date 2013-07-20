# Makefile

## Define paths to the libraries
LHAPDF=$(HOME)/work/Blackhole/lhapdf/lhapdf-5.8.9
PYTHIA8=$(HOME)/work/Blackhole/pythia8
INCLUDES=. /opt/local/include

## Compiler options
CC=g++
#CC=clang++
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

$(EXE): $(OBJS) tmp/PDFInterface.o src/main.cc
	$(CC) $(CCFLAGS) $(LDFLAGS) -o bin/$(EXE) src/main.cc $(OBJS) tmp/PDFInterface.o

clean:
	-rm -f bin/$(EXE)
	-rm -f tmp/*

tmp/PDFInterface.o: include/PDFInterface.h src/PDFInterface.cc
tmp/AbsModel.o: include/AbsModel.h src/AbsModel.cc
tmp/ADDModel.o: include/ADDModel.h src/ADDModel.cc
tmp/RSModel.o: include/RSModel.h src/RSModel.cc
tmp/Random.o: include/Random.h src/Random.cc
tmp/ConfigReader.o: include/ConfigReader.h src/ConfigReader.cc
tmp/Utility.o: include/Utility.h src/Utility.cc
# DO NOT DELETE

tmp/ADDModel.o: include/ADDModel.h include/AbsModel.h include/ConfigReader.h
tmp/ADDModel.o: include/PDFInterface.h include/Random.h
tmp/AbsModel.o: include/AbsModel.h include/ConfigReader.h
tmp/AbsModel.o: include/PDFInterface.h include/Random.h
tmp/ConfigReader.o: include/ConfigReader.h
tmp/PDFInterface.o: include/PDFInterface.h
tmp/RSModel.o: include/RSModel.h include/AbsModel.h include/ConfigReader.h
tmp/RSModel.o: include/PDFInterface.h include/Random.h
tmp/Random.o: include/Random.h /usr/include/boost/random.hpp
tmp/Utility.o: include/Utility.h
tmp/main.o: include/ConfigReader.h include/ADDModel.h include/AbsModel.h
tmp/main.o: include/PDFInterface.h include/Random.h
tmp/main.o: include/RSModel.h
tmp/main.o: include/Utility.h
