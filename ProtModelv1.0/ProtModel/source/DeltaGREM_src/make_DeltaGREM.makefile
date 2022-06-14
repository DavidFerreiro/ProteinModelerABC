#       Makefile
#
#       Copyright 2009 Unknown <fons@arnold.cbm.uam.es>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

####################

PROG=DeltaGREM


# Progam files

SRC=main_DeltaGREM.c REM.c random3.c Get_pars_DeltaGREM.c alignments.c \
	mutations.c Codes.c Input.c gen_code.c allocate.c output.c \
	read_pdb.c read.c Sec_str_all.c NeedlemanWunsch.c Profit_aux.c     

OBJ=main_DeltaGREM.o REM.o random3.o Get_pars_DeltaGREM.o alignments.o \
	mutations.o Codes.o Input.o gen_code.o allocate.o output.o \
	read_pdb.o read.o Sec_str_all.o NeedlemanWunsch.o Profit_aux.o  

# Compiler flags
CC=gcc
# CFLAGS= -O2
CFLAGS= -O2 -Wall -std=c99 -pedantic -g -pg -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
# -fbounds-check
LDFLAGS=-lm -lg2c
LDFLAGS=-lm -L/usr/lib64/libg2c.so.0.0.0

all: $(OBJ)
	$(CC) $(OBJ) -o $(PROG) $(LDFLAGS)

gnu: CC=gcc
gnu: LDFLAGS=-lm -lg2c
gnu: CFLAGS=-Wall -O3 -march=nocona -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
gnu: all

gnu-static: CC=gcc
gnu-static: LDFLAGS=-lm -lg2c -static
gnu-static: CFLAGS=-Wall -O3 -march=nocona -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE
gnu-static: all

intel: CC=icc
intel: LDFLAGS=-lm
intel: CFLAGS=-Wall -static -O3 -xHOST -ipo -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -vec-report
intel: all

intel-static: CC=icc
intel-static: LDFLAGS=-lm -static
intel-static: CFLAGS=-Wall -static -O3 -xHOST -ipo -D_FILE_OFFSET_BITS=64 -D_LARGE_FILE_SOURCE -vec-report
intel-static: all




%.o: %.cpp
	$(CC) $(LDFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm -fr $(OBJ) $(PROG)
