CXX = g++ -O0

CPLEX_HOME = /Users/sonndinh/Applications/IBM/ILOG/CPLEX_Studio1261/cplex
CONCERT_HOME = /Users/sonndinh/Applications/IBM/ILOG/CPLEX_Studio1261/concert
CPLEX_INCLUDE_PATH = $(CPLEX_HOME)/include
CONCERT_INCLUDE_PATH = $(CONCERT_HOME)/include
CPLEX_LIBS_PATH = $(CPLEX_HOME)/lib/x86-64_osx/static_pic
CONCERT_LIBS_PATH = $(CONCERT_HOME)/lib/x86-64_osx/static_pic

# Combined library path
CXXLN_LIBS_PATH = -L$(CPLEX_LIBS_PATH) -L$(CONCERT_LIBS_PATH)

# Compiler options
CXXOPT = -stdlib=libstdc++ -m64 -fPIC -fexceptions -DNDEBUG -DIL_STD

CXXFLAGS = $(CXXOPT) -I$(CPLEX_INCLUDE_PATH) -I$(CONCERT_INCLUDE_PATH)
CXXLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread

PROG = taskgen

all:: $(PROG) 

SRC = BlockingAnalysis.cpp \
	TaskSetGenerator.cpp \
	TestTaskSetGen.cpp

OBJS = BlockingAnalysis.o \
	TaskSetGenerator.o \
	TestTaskSetGen.o

$(OBJS): $(SRC)
	$(CXX) -c $(CXXFLAGS) $(SRC)

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(CXXLN_LIBS_PATH) $(CXXLNFLAGS) $(OBJS) -o $(PROG)

clean::
	rm -rf *.o
	rm taskgen
