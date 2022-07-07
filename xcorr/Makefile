# Makefile for fran

PROGRAM = xcorr
OBJECTS = xcorr.o

LIB = -lm -lgsl -lgslcblas -lsmurf -L${HOME}/Projects/lib_static
OPT = -I${HOME}/Projects/include -g #-funroll-all-loops -malign-double -fomit-frame-pointer -fexpensive-optimizations -Wno-deprecated 
FLAG = #-O3 -ffast-math -funroll-loops
CC = gcc
CXX = g++

${PROGRAM}: ${OBJECTS}
	${CXX} ${FLAG} ${OPT} ${OBJECTS} ${LIB} -o $(PROGRAM) 

install: ${OBJECTS}
	${CXX} ${FLAG} ${OPT} ${OBJECTS} ${LIB} -o $(PROGRAM) 
	cp -f ${PROGRAM} ${HOME}/Projects/bin

clean:
	rm -f ${PROGRAM} ${OBJECTS}

.cpp.o:
	$(CXX) $(OPT) -fPIC -o $@ -c $<

