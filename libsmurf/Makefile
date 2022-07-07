# Makefile for libsmurf

PROGRAM = libsmurf.so
A	= libsmurf.a
SO	= libsmurf.so

#JH remove xcorlib because it doesn't compile
OBJECTS = libsmurf.o fretdata.o hist.o traj.o xcorrlib.o langsim.o xtraj.o ltraj.o landscape.o friction.o ltraj_param.o
HEADERS = libsmurf.h fretdata.h hist.h traj.h xcorrlib.h langsim.h xtraj.h ltraj.h landscape.h friction.h dl_exception.h ltraj_param.h

#OBJECTS = libsmurf.o fretdata.o hist.o traj.o langsim.o xtraj.o ltraj.o landscape.o friction.o ltraj_param.o
#HEADERS = libsmurf.h fretdata.h hist.h traj.h langsim.h xtraj.h ltraj.h landscape.h friction.h dl_exception.h ltraj_param.h


CROOTPATH = /usr/share/CLAPACK
INCDIRS = -I$(CROOTPATH)/SRC -I$(CROOTPATH) -I$(CROOTPATH)/INCLUDE
F2CDIR  = $(CROOTPATH)/F2CLIBS
LDLIBS  = $(CROOTPATH)/lapack_LINUX.a \
          $(CROOTPATH)/blas_LINUX.a \
	$(F2CDIR)/libF77.a $(F2CDIR)/libI77.a -lm -lgsl -lgslcblas
LIB = $(LDLIBS) -L${HOME}/Projects/lib/ -lm -lgsl -lgslcblas
OPT = $(INCDIRS) -g -pg -fPIC -Wall # -funroll-all-loops -malign-double -fomit-frame-pointer -fexpensive-optimizations
CXXFLAGS = #-O3 -ffast-math -funroll-loops -Wno-deprecated
CC = gcc $(INCDIRS)
CXX= g++
AR= ar rcs

so: $(OBJECTS) ${HEADERS}
	${CXX} ${OPT} ${CXXFLAGS} -shared -o ${PROGRAM} ${OBJECTS} ${LIB}

a: ${OBJECTS} ${HEADERS}
	${AR} ${A} ${OBJECTS}

install: so a
#	rm -f ${HOME}/Projects/lib/libsmurf.a
	cp libsmurf.so ${HOME}/Projects/lib/
	cp libsmurf.a ${HOME}/Projects/lib_static/
	cp ${HEADERS} ${HOME}/Projects/include/

installa: a
	rm -f ${HOME}/Projects/lib/libsmurf.so
	cp libsmurf.a ${HOME}/Projects/lib_static/
	cp ${HEADERS} ${HOME}/Projects/include/


test: test.cpp
	g++ test.cpp -o test -lm -lsmurf ${CXXFLAGS}

.cpp.o:
	$(CXX) $(OPT) -fPIC -o $@ -c $<

clean:
	rm -f ${A} ${SO} ${OBJECTS}
