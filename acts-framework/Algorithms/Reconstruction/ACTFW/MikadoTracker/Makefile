USEROOT = yes
ISROOT := $(shell which root-config)

ifndef ISROOT
  USEROOT = no
endif

CC=g++
CFLAGS=-c -Wall -std=c++14 -O3 -pthread
INCLUDES=-I. -I./analysis 

ifeq ($(USEROOT), yes)
INCLUDES:=$(INCLUDES) -I$(shell root-config --incdir)
CFLAGS:=$(CFLAGS) -DUSEROOT
LDFLAGS:=$(shell root-config --glibs)
endif

#CFLAGS:=$(CFLAGS) -DSGSafeRun

#CFLAGS:=$(CFLAGS) -DSGUseMCInReco

#(root-config --incdir)
LIBDIRS=-L/usr/lib
LIBS=
LDFLAGS:=$(LDFLAGS) -pthread

SRC=TrackModelPhysical.cxx SearchLayer.cxx  Tracker.cxx reconstruction.cxx Geo.cxx Cuts.cxx \
analysis/PolynomFit.cxx DataStructures.cxx Learning.cxx AccuracyEvaluator.cxx Engine.cxx \
TrackSelector.cxx EventReader.cxx

SOURCES=$(SRC) analysis/TrackFitTest.cxx  analysis/AnalyseGeometry.cxx analysis/AnalyseField.cxx 

HEADERS=$(SRC:.cxx=.h) util.h
OBJECTS=$(SOURCES:.cxx=.o)
EXECUTABLE=reco

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) $(LIBS) -o $@

%.o: %.cxx $(HEADERS) Makefile
	$(CC) $(CFLAGS) $(INCLUDES) $< -o $@ 


clean:
	rm -f reco $(OBJECTS) *~ mysubmission.csv \#*; 
#	cd phase2; . ./clean.sh 

cleanall:
	rm reco $(OBJECTS) *.~ \#* *.root
