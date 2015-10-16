SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic

ifeq ("$(HOSTNAME)","phdstulin2.isye.gatech.edu")
	HPCOPTPATH = /hpcopt
else
	HPCOPTPATH = /opt
endif

CPLEXDIR      = $(HPCOPTPATH)/ibm/ILOG/CPLEX_Studio126/CPLEX_Studio/cplex
CONCERTDIR    = $(HPCOPTPATH)/ibm/ILOG/CPLEX_Studio126/CPLEX_Studio/concert

CCC = g++ -std=c++0x

CCOPT = -O3 -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD

CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)

CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread

CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 

OBJ_DIR = obj
OUT_RELEASE = bin/gep

OBJ_FILE = $(OBJ_DIR)/mt64.o $(OBJ_DIR)/functions.o $(OBJ_DIR)/siddp.o

all: siddp

siddp: $(OBJ_FILE)
	$(CCC) $(CCFLAGS) $(OBJ_FILE) -o $(OUT_RELEASE) $(CCLNFLAGS)

$(OBJ_DIR)/siddp.o: siddp.cpp
	$(CCC) -c $(CCFLAGS) siddp.cpp -o $(OBJ_DIR)/siddp.o

$(OBJ_DIR)/mt64.o: mt64.cpp
	$(CCC) -c $(CCFLAGS) mt64.cpp -o $(OBJ_DIR)/mt64.o

$(OBJ_DIR)/functions.o: functions.cpp
	$(CCC) -c $(CCFLAGS) functions.cpp -o $(OBJ_DIR)/functions.o

