
all: mmpi

ROOT_INCLUDES = -I $(shell root-config --incdir)
ROOT_LIBS =  -L/usr/lib64/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix  -lMathCore -lThread -pthread -lm -ldl -rdynamic

mmpi: mmpi.cc murmur.h
	mpicxx -std=c++11 -Wall -O2 $(ROOT_INCLUDES) $(ROOT_LIBS) -o mmpi mmpi.cc 

