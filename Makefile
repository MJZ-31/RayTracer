CC=g++
CFLAGS=-std=c++20
CLIBS=-lpng

SRCDIR=src
BUILDDIR=build
OBJDIR=$(BUILDDIR)/obj
BINDIR=$(BUILDDIR)/bin

all:
	$(CC) $(CFLAGS) src/main.cpp $(CLIBS)