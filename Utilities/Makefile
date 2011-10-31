# Title: Makefile for Python library files KPDF.c, Cmap.c and Cstats.c
# Author: Craig Arthur, craig.arthur@ga.gov.au
# CreatedOn: 2009-12-18
# Description: Makefile for compiling C versions of a number of Python modules
#              This version is for use in *nix operating systems that have
#              GNU make installed.
#
# Version: $Rev$
# 
# $Id$

CC 	= gcc
VPATH	= /usr/local/python/include/python2.5:/usr/local/python/lib/python2.5
CFLAGS 	= -Wall -O -fPIC -std=c99
CFLAGS	+= $(patsubst %,-I%,$(subst :, ,$(VPATH)))
OBJS	= KPDF.o Cstats.o Cmap.o
LDFLAGS = -shared

all: $(OBJS)

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@
	$(CC) $(LDFLAGS) $@ -o $(@:.o=.so)

clean: 
	rm $(OBJS) $(OBJS:.o=.so)
