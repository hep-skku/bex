# Makefile

## Compiler options
## Detect OSX
ifeq ($(shell uname),Darwin)
CC=clang++
#CC=g++
SWPATH=$(HOME)/sw/osx_x86_64_gcc482
BOOSTFLAGS=-L/opt/local/lib -lboost_iostreams-mt
else
CC=g++
SWPATH=$(HOME)/sw/slc5_x86_64_gcc412
BOOSTFLAGS=-lboost_iostreams
endif
INCLUDES=. /opt/local/include

## Define paths to the libraries
LHAPDF=$(SWPATH)/lhapdf
PYTHIA8=$(SWPATH)/pythia8
HEPMC=$(SWPATH)/HepMC
MYPATH=$(shell pwd)

EXE=bex
CCFLAGS=-Wall $(addprefix -I,$(INCLUDES)) -O3 -DBASEDIR=\"$(shell pwd)\"
LDFLAGS=-L$(LHAPDF)/lib -lLHAPDF -lm -lz $(BOOSTFLAGS)

## Detect ROOT for debugging
ifdef ROOTSYS
CCFLAGS:=$(CCFLAGS) $(shell root-config --cflags) -DDEBUGROOT
LDFLAGS:=$(LDFLAGS) $(shell root-config --ldflags --libs) -DDEBUGROOT
endif

## Actions
all: $(EXE) SCRIPTS

EXCLUDES=main.cc PDFInterface.cc pythiaHadronizer.cc
SRCS=$(filter-out $(EXCLUDES),$(notdir $(wildcard src/*.cc)))
OBJS=$(addprefix tmp/,$(SRCS:.cc=.o))

$(OBJS):
	$(CC) $(CCFLAGS) -o $@ -c src/$(*F).cc

tmp/PDFInterface.o:
	$(CC) $(CCFLAGS) -I$(LHAPDF)/include -o $@ -c src/$(*F).cc

$(EXE): $(OBJS) tmp/PDFInterface.o src/main.cc
	$(CC) $(CCFLAGS) $(LDFLAGS) -o bin/$(EXE) src/main.cc $(OBJS) tmp/PDFInterface.o

pythia: src/pythiaHadronizer.cc
	$(CC) src/pythiaHadronizer.cc -o bin/hadronizer $(CCFLAGS) $(LDFLAGS) \
        -I$(PYTHIA8)/include -L$(PYTHIA8)/lib/archive -lpythia8 -lhepmcinterface \
        -I$(HEPMC)/include -L$(HEPMC)/lib -lHepMC -L$(HEPMC)/lib

SCRIPTS:
	@# bex run script
	@echo "#!/bin/bash" > bin/run.sh
	@echo 'export LD_LIBRARY_PATH=$${LD_LIBRARY_PATH}:'$(LHAPDF)/lib >> bin/run.sh
	@echo $(MYPATH)/bin/$(EXE) '$$@' >> bin/run.sh
	@chmod +x bin/run.sh
	@# pythiaHadronizer run script
	@echo "#!/bin/bash" > bin/runHadronizer.sh
	@echo 'export LD_LIBRARY_PATH=$${LD_LIBRARY_PATH}:'$(LHAPDF)/lib:$(HEPMC)/lib:$(PYTHIA8)/lib >> bin/runHadronizer.sh
	@echo export PYTHIA8DATA=${PYTHIA8}/xmldoc >> bin/runHadronizer.sh
	@echo $(MYPATH)/bin/hadronizer '$$@' >> bin/runHadronizer.sh
	@chmod +x bin/runHadronizer.sh

clean:
	-rm -f bin/*
	-rm -f tmp/*

tmp/PDFInterface.o: include/PDFInterface.h src/PDFInterface.cc
tmp/AbsModel.o: include/AbsModel.h src/AbsModel.cc
tmp/ADDModel.o: include/ADDModel.h src/ADDModel.cc
tmp/RSModel.o: include/RSModel.h src/RSModel.cc
tmp/Random.o: include/Random.h src/Random.cc
tmp/ConfigReader.o: include/ConfigReader.h src/ConfigReader.cc
tmp/Utility.o: include/Utility.h src/Utility.cc
tmp/NVector.o: include/NVector.h src/NVector.cc
# DO NOT DELETE

tmp/ADDModel.o: include/AbsModel.h include/ConfigReader.h
tmp/ADDModel.o: include/PDFInterface.h include/Random.h
tmp/AbsModel.o: include/ConfigReader.h
tmp/AbsModel.o: include/PDFInterface.h include/Random.h
tmp/RSModel.o: include/AbsModel.h include/ConfigReader.h
tmp/RSModel.o: include/PDFInterface.h include/Random.h
tmp/main.o: include/ConfigReader.h include/ADDModel.h include/AbsModel.h
tmp/main.o: include/PDFInterface.h include/Random.h
tmp/main.o: include/RSModel.h
tmp/main.o: include/Utility.h
