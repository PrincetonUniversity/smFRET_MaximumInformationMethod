PROGRAM1 = boots
OBJECTS1 = boots.o
HEADERS = boots.h

CROOTPATH = /usr/share/CLAPACK
INCDIRS = -I$(CROOTPATH)/SRC -I$(CROOTPATH)
F2CDIR = $(CROOTPATH)/F2CLIBS
LDLIBS = -l$(CROOTPATH)/lapack_LINUX.a \
	 -l$(CROOTPATH)/blas_LINUX.a \
	 -l$(CROOTPATH)/libF77.a -l$(F2CDIR)/libI77.a
LIB = $(LBLIBS) -lm -lgsl -lgslcblas  -lsmurf -L../libsmurf/ #-lsmurf -L${HOME}/Projects/lib/ -L${HOME}/Projects/lib_static/
OPT = $(INCDIRS) -I../libsmurf #-I${HOME}/Projects/include -g -pg# -funroll-all-loops -malign-double -fomit-frame-pointer -fexpensive-optimizations -Wno-deprecated 
FLAG = #-O3 -ffast-math -funroll-loops
CC = gcc $(INCDIRS)
CXX = g++

${PROGRAM1}: ${OBJECTS1} ${HEADERS}
	${CXX} ${OPT} ${FLAG} -o ${PROGRAM1} ${OBJECTS1} ${LIB} 

install: ${OBJECTS1} ${HEADERS}
	${CXX} ${OPT} ${FLAG} -o ${PROGRAM1} ${OBJECTS1} ${LIB} 
	cp -f ${PROGRAM1} ${HOME}/Projects/bin

clean:
	rm -f ${PROGRAM1} ${OBJECTS1}

.cpp.o: $<
	$(CXX) $(OPT) -fPIC -o $@ -c $<
