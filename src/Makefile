#****************************************************
# (c) Copyright 2001, All rights reserved.
#
# gmhmme3 - GeneMark.hmm eukaryotic version 3
#   Alex Lomsadze and Mark Borodovsky
#
# Makefile for GeneMark.hmm eukaryotic
#****************************************************

srcdir	= .
CXX	= g++
LDFLAGS	= -lstdc++ -static
LIBS	= -lm
CFLAGS	= -O3

SRCS	= main.cpp cmdline.cpp astring.cpp fload.cpp sequence.cpp \
	map.cpp  markov.cpp matrix.cpp model.cpp \
	parser.cpp path.cpp precalc.cpp site.cpp viterbi.cpp bp.cpp

OBJS	= main.o cmdline.o astring.o fload.o sequence.o \
	map.o markov.o matrix.o model.o \
	parser.o path.o precalc.o site.o viterbi.o bp.o
 	  
HDRS	= common.h cmdline.h astring.h fload.h alist.h sequence.h \
	map.h baseviterbi.h length.h markov.h matrix.h model.h \
	parser.h path.h precalc.h site.h viterbi.h bp.h

.SUFFIXES:	.o .cpp ${SUFFIXES}

.cpp.o:
	${CXX} -o $@ -c $< -I${srcdir} ${CFLAGS}

all:	gmhmme3

gmhmme3:	${OBJS} ${HDRS}
	${CXX} -o $@ ${OBJS} ${LDFLAGS} ${LIBS}
	strip gmhmme3

foretp:	
	${MAKE} gmhmme3 CFLAGS="${CFLAGS} -DETP__"
	strip gmhmme3

clean:
	rm -f ${OBJS}
	rm -f gmhmme3

