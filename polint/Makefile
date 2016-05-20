###############################################################################
#   This Makefile shows how to compile all C++, C and Fortran
#   files found in $(SRCDIR) directory.
#   Linking is done with g++. Need to have $ROOTSYS defined
###############################################################################

###############################################################################
MYOS        := $(shell uname)
ARCH        := $(shell uname -m)
SHELL       := /bin/bash

###############################################################################
EXEFILE1    := polin
EXEFILE2    := polel

SRCDIR      := src
INCDIR      := src
OBJDIR      := obj.$(ARCH)

MAINFILE1   := polin.f
MAINFILE2   := polel.f

INCDIRS     := $(INCDIR)
MAINOBJ1    := $(OBJDIR)/$(addsuffix .o,$(basename $(MAINFILE1)))
MAINOBJ2    := $(OBJDIR)/$(addsuffix .o,$(basename $(MAINFILE2)))

###############################################################################
# Compiler
AR          := ar
CC          := gcc
CXX         := g++
FC          := gfortran
LD          := g++

###############################################################################
# Flags
ifeq ($(ARCH),i686)
    MODE    := -m32
else
    MODE    := -m64
endif
INCFLAGS    := $(addprefix -I,$(INCDIRS))
CFLAGS      := -Wall -fPIC -Wno-unused-variable -g $(MODE) $(INCFLAGS)
CXXFLAGS    := -Wall -fPIC -Wno-unused-variable -g $(MODE) $(INCFLAGS)
FFLAGS      := -Wall -fPIC -Wno-unused-variable -g $(MODE) $(INCFLAGS)
ifeq ($(MYOS),Darwin) #in Darwin, do not use -fno-leading-underscore
    FFLAGS  += -fno-second-underscore -fbounds-check \
               -ffixed-line-length-none
else
    FFLAGS  += -fno-leading-underscore -fno-second-underscore \
               -fbounds-check -ffixed-line-length-none
endif
CPPFLAGS    :=
LDFLAGS     := -g $(MODE)

###############################################################################
# Generate file lists
# source files
FSOURCES    := $(wildcard $(SRCDIR)/*.[Ff])
CSOURCES    := $(wildcard $(SRCDIR)/*.[Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][Cc])
CSOURCES    += $(wildcard $(SRCDIR)/*.[Cc][XxPp][XxPp])
SOURCES     := $(FSOURCES) $(CSOURCES)

# header files
HEADERS     := $(wildcard $(INCDIR)/*.[Hh])
HEADERS     += $(wildcard $(INCDIR)/*.[Hh][Hh])

# obj and dep files
OBJS        := $(addsuffix .o,$(basename $(SOURCES)))
OBJS        := $(patsubst $(SRCDIR)/%.o,$(OBJDIR)/%.o,$(OBJS))
DEPS        := $(subst .o,.d,$(OBJS))

###############################################################################
# Libraries
SYSLIBS     := -lgfortran
LIBS        += $(SYSLIBS)

###############################################################################
# Pre-defined commands
define make-dependency-cc
echo "Making dependency for file $< ......"
set -e; [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR); \
$(CC) -MM $(CPPFLAGS) $(CFLAGS) $< > $@.$$$$; \
sed 's,\($*\)\.o[ :]*,$(OBJDIR)/\1.o $@ : ,g' < $@.$$$$ > $@; \
rm -f $@.$$$$
endef

define make-dependency-cxx
echo "Making dependency for file $< ......"
set -e; [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR); \
$(CXX) -MM $(CPPFLAGS) $(CXXFLAGS) $< > $@.$$$$; \
sed 's,\($*\)\.o[ :]*,$(OBJDIR)/\1.o $@ : ,g' < $@.$$$$ > $@; \
rm -f $@.$$$$
endef

define compile-cc
echo "Compiling $< ......"
$(CC) $(CFLAGS) -c $< -o $@
endef

define compile-cxx
echo "Compiling $< ......"
$(CXX) $(CXXFLAGS) -c $< -o $@
endef

define compile-fc
echo "Compiling $< ......"
set -e; [ -d $(OBJDIR) ] || mkdir -p $(OBJDIR); \
$(FC) $(FFLAGS) -c $< -o $@
endef

###############################################################################
# Targets
VPATH       := $(SRCDIR)

.SUFFIXES:
.SUFFIXES: .c .C .cc .CC .cpp .cxx .f .F
.SUFFIXES: .h .hh .H .HH
.SUFFIXES: .o .d

.PHONY: all exe1 exe2
.PHONY: clean distclean
.PHONY: test help env

###############################################################################
all: exe1 exe2

###############################################################################
exe1: $(OBJS) $(MAINOBJ1)
	@$(LD) $(LDFLAGS) -o $(EXEFILE1) $(OBJS) $(MAINOBJ1) $(LIBS)
	@echo "Linking $(EXEFILE1) ...... done!"

exe2: $(OBJS) $(MAINOBJ2)
	@$(LD) $(LDFLAGS) -o $(EXEFILE2) $(OBJS) $(MAINOBJ2) $(LIBS)
	@echo "Linking $(EXEFILE2) ...... done!"

$(MAINOBJ1): $(MAINFILE1)
	@$(compile-fc)

$(MAINOBJ2): $(MAINFILE2)
	@$(compile-fc)

###############################################################################
# Make dependencies
$(OBJDIR)/%.d: %.c
	@$(make-dependency-cc)

$(OBJDIR)/%.d: %.C
	@$(make-dependency-cxx)

$(OBJDIR)/%.d: %.cc
	@$(make-dependency-cxx)

$(OBJDIR)/%.d: %.CC
	@$(make-dependency-cxx)

$(OBJDIR)/%.d: %.cpp
	@$(make-dependency-cxx)

$(OBJDIR)/%.d: %.cxx
	@$(make-dependency-cxx)

ifneq ($(DEPS),)
    sinclude $(DEPS)
endif

###############################################################################
# Compile source files
$(OBJDIR)/%.o: %.c
	@$(compile-cc)

$(OBJDIR)/%.o: %.C
	@$(compile-cxx)

$(OBJDIR)/%.o: %.cc
	@$(compile-cxx)

$(OBJDIR)/%.o: %.CC
	@$(compile-cxx)

$(OBJDIR)/%.o: %.cpp
	@$(compile-cxx)

$(OBJDIR)/%.o: %.cxx
	@$(compile-cxx)

$(OBJDIR)/%.o: %.f
	@$(compile-fc)

$(OBJDIR)/%.o: %.F
	@$(compile-fc)

###############################################################################
# Clean
clean:
	@rm -f $(EXEFILE1) $(EXEFILE2)
	@rm -f $(OBJDIR)/*
	@rm -f *~ *# */*~ */*#

distclean: clean

test:
	@echo "MYOS:     $(ARCH) $(MYOS)"
	@echo "CFLAGS:   $(CFLAGS)"
	@echo "CXXFLAGS: $(CXXFLAGS)"
	@echo "FFLAGS:   $(FFLAGS)"
	@echo "LDFLAGS:  $(LDFLAGS)"
	@echo "CSOURCES: $(CSOURCES)"
	@echo "FSOURCES: $(FSOURCES)"
	@echo "HEADERS:  $(HEADERS)"
	@echo "OBJS:     $(OBJS)"
	@echo "DEPS:     $(DEPS)"

help: test

env: test
