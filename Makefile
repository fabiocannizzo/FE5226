# Possible arguments
# DEBUG=1               [ 0 | 1 ],         default 0
# COMPILER              [ g++ | clang++ ], default: g++ on Linux
#

DEBUG ?= 0

# platform
PLATFORM := $(shell uname -s)
$(info PLATFORM: $(PLATFORM))

# compiler
COMPILER ?= g++
CC = $(COMPILER)
$(info COMPILER: $(CC))

# directories
BINDIR := bin
ifeq ($(DEBUG),1)
  BINDIR := $(BINDIR)-dbg
endif
$(info BINDIR: $(BINDIR))

# file extensions
ifeq ($(PLATFORM), Linux)
  OBJ=.o
  EXE=.out
  DLL=.so
  LIB=.a
else
# cygwin or mingw
  OBJ=.obj
  EXE=.exe
  DLL=.dll
  LIB=.lib
endif
$(info EXTENSIONS:  $(OBJ) $(EXE) $(DLL) $(LIB))

# detect cpp files containing main and filter for files tracked in git
GIT_FILES := $(shell git ls-files 2>/dev/null)
MAINS := $(shell grep -l "int main" *.cpp)
ifneq ($(GIT_FILES),)
  MAINS := $(filter $(GIT_FILES),$(MAINS))
endif
$(info MAINS:   $(MAINS))

# target variable containing all compiled filenames
TARGETS := $(patsubst %.cpp,$(BINDIR)/%$(OBJ),$(MAINS))
$(info TARGETS: $(TARGETS))

DEPFILES := $(patsubst %.cpp,$(BINDIR)/%.d,$(MAINS))

DEPFLAGS = -MT $@ -MMD -MP -MF $(BINDIR)/$*.d
CFLAGS := -c -std=c++20 -march=native -Wall

ifeq ($(DEBUG),1)
   CFLAGS += -g -DDEBUG
else
   CFLAGS += -O3
endif

.PHONY: all clean MC

all: $(TARGETS)
	$(MAKE) -C MC

MC:
	$(MAKE) -C MC

.SECONDARY: $(TARGETS) $(DEPFILES)

$(BINDIR)/%$(OBJ): %.cpp Makefile $(BINDIR)/%.d | $(BINDIR)
	$(CC) $(CFLAGS) $(DEPFLAGS) $< -o $@

$(DEPFILES):

include $(wildcard $(DEPFILES))

# clean
clean:
	rm -rf $(BINDIR) $(BINDIR)-dbg
	$(MAKE) -C MC clean

# create bin folder
$(BINDIR):
	mkdir -p $@
