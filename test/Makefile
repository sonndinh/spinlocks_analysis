CXX_HOME = /export/austen/home/son/tools/gcc-4.9.2-install
CXX = $(CXX_HOME)/bin/g++ -O0

CPLEX_HOME = /export/austen/home/son/tools/ibm/ILOG/CPLEX_Studio1261/cplex
CONCERT_HOME = /export/austen/home/son/tools/ibm/ILOG/CPLEX_Studio1261/concert
CPLEX_INCLUDE_PATH = $(CPLEX_HOME)/include
CONCERT_INCLUDE_PATH = $(CONCERT_HOME)/include
CPLEX_LIBS_PATH = $(CPLEX_HOME)/lib/x86-64_linux/static_pic
CONCERT_LIBS_PATH = $(CONCERT_HOME)/lib/x86-64_linux/static_pic

# Combined library path
CXXLN_LIBS_PATH = -L$(CPLEX_LIBS_PATH) -L$(CONCERT_LIBS_PATH)

# Compiler options
CXXOPT = -std=c++11 -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD

CXXFLAGS = $(CXXOPT) -I$(CPLEX_INCLUDE_PATH) -I$(CONCERT_INCLUDE_PATH)
CXXLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread

PROG = test_prng test_miqcp

all:: $(PROG)

SRC = TestRanNumGen.cpp \
	TestMIQCP.cpp

OBJS = TestRanNumGen.o \
	TestMIQCP.o

$(OBJS): $(SRC)
	$(CXX) $(CXXFLAGS) -c $(SRC)

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(CXXLN_LIBS_PATH) TestRanNumGen.o -o test_prng $(CXXLNFLAGS)
	$(CXX) $(CXXFLAGS) $(CXXLN_LIBS_PATH) TestMIQCP.o -o test_miqcp $(CXXLNFLAGS)

clean::
	rm -rf *.o
	rm test_prng test_miqcp
