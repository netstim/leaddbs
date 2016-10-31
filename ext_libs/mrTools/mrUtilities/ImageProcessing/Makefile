MLAB = /Applications/MATLAB7.7.app/

MXSFX = mexmaci64
MEX = ${MLAB}/bin/mex

MFLAGS = -largeArrayDims -cxx -maci64
INC = -I ${MLAB}/extern/include
LIB = -L ${MLAB}/extern/lib

CC = cc -Wall -pedantic -no-cpp-precomp -arch x86_64 
C_OPTIMIZE_SWITCH = -O2    ## For GCC
CFLAGS = ${C_OPTIMIZE_SWITCH} ${INC} ${LIB}

all: corrDn.${MXSFX} upConv.${MXSFX} 

clean:
	/bin/rm *.o

corrDn.${MXSFX}: corrDn.o wrap.o convolve.o edges.o
	${MEX} ${MFLAGS} corrDn.o wrap.o convolve.o edges.o

upConv.${MXSFX}: upConv.o wrap.o convolve.o edges.o
	${MEX} ${MFLAGS} upConv.o wrap.o convolve.o edges.o

convolve.o wrap.o edges.o: convolve.h 

%.o : %.c
	${CC} -c ${CFLAGS} $<		

