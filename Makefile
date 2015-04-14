# Makefile is a part of the PYTHIA event generator.
# Copyright (C) 2015 Torbjorn Sjostrand.
# PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
# Please respect the MCnet Guidelines, see GUIDELINES for details.
# Author: Philip Ilten, September 2014.
#
# This is is the Makefile used to ILC hz simulation (based on PYTHIA examples 
# Makefile )
# Example usage is:
#     make ilchz
# For help using the make command please consult the local system documentation,
# i.e. "man make" or "make --help".

################################################################################
# VARIABLES: Definition of the relevant variables from the configuration script.
################################################################################

# Include the configuration.
PREFIX_BIN=/opt/pythia8205/bin
PREFIX_INCLUDE=/opt/pythia8205/include
PREFIX_LIB=/opt/pythia8205/lib
PREFIX_SHARE=/opt/pythia8205/share/Pythia8
# Compilation flags (see ./configure --help for further documentation).
ENABLE_SHARED=true
CXX=g++
CXX_COMMON=-O2 -ansi -pedantic -W -Wall -Wshadow -fPIC
CXX_SHARED=-shared
CXX_SONAME=-Wl,-soname
LIB_SUFFIX=.so

# HEPMC2 configuration.
HEPMC2_USE=true
HEPMC2_BIN=/usr/bin/
HEPMC2_INCLUDE=/usr/include
HEPMC2_LIB=/usr/lib

# ROOT configuration.
ROOT_USE=true
ROOT_BIN=/usr/bin/
ROOT_INCLUDE=/usr/include/root
ROOT_LIB=/usr/lib

# Check distribution (use local version first, then installed version).
ifneq ("$(wildcard ../lib/libpythia8.a)","")
  PREFIX_LIB=../lib
  PREFIX_INCLUDE=../include
endif
CXX_COMMON:=-I$(PREFIX_INCLUDE) $(CXX_COMMON) -Wl,-rpath $(PREFIX_LIB) -ldl

################################################################################
# RULES: Definition of the rules used to build the PYTHIA examples.
################################################################################

# Rules without physical targets (secondary expansion for specific rules).
.SECONDEXPANSION:
.PHONY: all clean

# All targets (no default behavior).
all:
	@echo "Usage: make ilchz"

# The Makefile configuration.
# PYTHIA libraries.
$(PREFIX_LIB)/libpythia8.a :
	$(error Error: PYTHIA must be built, please run "make"\
                in the top PYTHIA directory)

ilchz: $$@.cc $$@.h 
	@echo "Compiling ilchz executable..."
	@export LD_LIBRARY_PATH=$$LD_LIBRARY_PATH:$(ROOT_LIB);\
	$(CXX) $^ -o $@ -w -I$(ROOT_INCLUDE) -I$(HEPMC2_INCLUDE) $(CXX_COMMON)\
	 -Wl,-rpath $(ROOT_LIB) -L$(PREFIX_LIB) -L$(HEPMC2_LIB) `$(ROOT_BIN)root-config --glibs` \
	 -lHepMC -lpythia8

# Clean.
clean:
	@rm -f ilchz;
