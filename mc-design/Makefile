# Every file exNN-<Name>.cpp is a stand-alone example: one file, one program.
#
# make                 build all examples          -> bin/exNN-<Name>.out
# make ex04-EuroCallMarketVec   build one example
# make run             build and run all examples
# make clean
#
# Possible arguments
# DEBUG=1               [ 0 | 1 ],         default 0
# COMPILER              [ g++ | clang++ ], default: g++

DEBUG ?= 0
COMPILER ?= g++
CC = $(COMPILER)

PLATFORM := $(shell uname -s)
ifeq ($(PLATFORM), Linux)
  EXE=.out
else
  EXE=.exe
endif

BINDIR := bin
ifeq ($(DEBUG),1)
  BINDIR := $(BINDIR)-dbg
endif

CFLAGS := -std=c++20 -march=native -Wall -Werror
ifeq ($(DEBUG),1)
  CFLAGS += -g -DDEBUG
else
  CFLAGS += -O3
endif

EXAMPLES := $(sort $(patsubst %.cpp,%,$(wildcard ex*.cpp)))
TARGETS  := $(patsubst %,$(BINDIR)/%$(EXE),$(EXAMPLES))

all : $(TARGETS)

# one program per source file
$(BINDIR)/%$(EXE) : %.cpp Makefile | $(BINDIR)
	$(CC) $(CFLAGS) $< -o $@

# make ex04-EuroCallMarketVec
$(EXAMPLES) : % : $(BINDIR)/%$(EXE)
.PHONY: $(EXAMPLES)

.PHONY: run
run : $(TARGETS)
	@for t in $(TARGETS); do echo "=== $$t"; ./$$t || { echo "FAILED: $$t"; exit 1; }; done

.PHONY: clean
clean :
	rm -rf bin bin-dbg

$(BINDIR) :
	mkdir -p $@
