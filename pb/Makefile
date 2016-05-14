########################################################################
#   This Makefile shows how to compile all C++, C and Fortran
#   files found in $(SRCDIR) directory.
#   Linking is done with g++. Need to have $ROOTSYS defined
########################################################################

########################################################################
MYOS        := $(shell uname)
ARCH        := $(shell uname -m)
SHELL       := /bin/bash

########################################################################
EXEFILE     := pbF12

########################################################################
# File list
OBJS        := F1F209.o F1F209Wrapper.o pbF12.o

########################################################################
# Compiler
AR          := ar
CC          := gcc
CXX         := g++
FF          := gfortran
LD          := g++

########################################################################
# Flags
ifeq ($(ARCH),i686)
    MODE    := -m32
else
    MODE    := -m64
endif
CFLAGS      := -Wall -fPIC -Wno-unused-variable -O3 -g $(MODE)
CXXFLAGS    := -Wall -fPIC -Wno-unused-variable -O3 -g $(MODE)
FFLAGS      := -Wall -fPIC -Wno-unused-variable -O3 -g $(MODE)
ifeq ($(MYOS),Darwin) #in Darwin, do not use -fno-leading-underscore
    FFLAGS  += -fno-second-underscore -fno-automatic -fbounds-check \
               -funroll-all-loops -fdollar-ok -ffixed-line-length-none \
               -fno-range-check
else
    FFLAGS  += -fno-leading-underscore -fno-second-underscore \
               -fno-automatic -fbounds-check -funroll-all-loops \
               -fdollar-ok -ffixed-line-length-none -fno-range-check
endif
CPPFLAGS    :=
LDFLAGS     := -O3 -g $(MODE)

########################################################################
# Libs
SYSLIBS     := -lstdc++ -lgfortran

########################################################################
# You can specify the .SUFFIXES
.SUFFIXES: .c .C .cc .CC .cpp .cxx .f .F
.PHONY: all clean test

########################################################################
all: exe

########################################################################
exe: $(OBJS)
	@$(LD) $(LDFLAGS) -o $(EXEFILE) $(OBJS) $(SYSLIBS)
	@echo "Linking $(EXEFILE) ...... done!"

########################################################################
%.o: %.c
	@echo Compiling $< ......
	@$(CC) -c $< -o $@ $(CFLAGS)

%.o: %.C
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@ $(CXXFLAGS)

%.o: %.cc
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@ $(CXXFLAGS)

%.o: %.CC
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@ $(CXXFLAGS)

%.o: %.cpp
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@ $(CXXFLAGS)

%.o: %.cxx
	@echo Compiling $< ......
	@$(CXX) -c $< -o $@ $(CXXFLAGS)

%.o: %.f
	@echo Compiling $< ......
	@$(FF) -c $< -o $@ $(FFLAGS)

%.o: %.F
	@echo Compiling $< ......
	@$(FF) -c $< -o $@ $(FFLAGS)

########################################################################
clean:
	@rm -f *.o
	@rm -f *~ *# */*~ */*#
	@rm -f $(EXEFILE)

test:
	@echo "MYOS:     $(ARCH) $(MYOS)"
	@echo "CFLAGS:   $(CFLAGS)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "FFLAGS:   $(FFLAGS)"
	@echo "LDFLAGS:  $(LDFLAGS)"
	@echo "OBJS:     $(OBJS)"

help: test

env: test
