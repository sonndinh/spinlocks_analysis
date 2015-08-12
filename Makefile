CXX_HOME = /export/austen/home/son/tools/gcc-4.9.2-install
CXX = $(CXX_HOME)/bin/g++ -O0

### Cplex paths
#CPLEX_HOME = /Users/sonndinh/Applications/IBM/ILOG/CPLEX_Studio1261/cplex
#CONCERT_HOME = /Users/sonndinh/Applications/IBM/ILOG/CPLEX_Studio1261/concert
CPLEX_HOME = /export/austen/home/son/tools/ibm/ILOG/CPLEX_Studio1261/cplex
CONCERT_HOME = /export/austen/home/son/tools/ibm/ILOG/CPLEX_Studio1261/concert
CPLEX_INCLUDE_PATH = $(CPLEX_HOME)/include
CONCERT_INCLUDE_PATH = $(CONCERT_HOME)/include
CPLEX_LIBS_PATH = $(CPLEX_HOME)/lib/x86-64_linux/static_pic
CONCERT_LIBS_PATH = $(CONCERT_HOME)/lib/x86-64_linux/static_pic

### SCIP paths
SCIP_HOME = /export/austen/home/son/tools/scip-3.1.1-install
SCIP_INCLUDE_PATH = $(SCIP_HOME)/include
SCIP_LIBS_PATH = $(SCIP_HOME)/lib

### ZIMPL paths
ZIMPL_HOME = /export/austen/home/son/tools/scipoptsuite-3.1.1/zimpl-3.3.2
ZIMPL_LIBS_PATH = $(ZIMPL_HOME)/lib

# Combined library path
CXXLN_LIBS_PATH = -L$(CPLEX_LIBS_PATH) -L$(CONCERT_LIBS_PATH) -L$(SCIP_LIBS_PATH) -L$(ZIMPL_LIBS_PATH)

# Compiler options
CXXOPT = -std=c++11 -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD

CXXFLAGS = $(CXXOPT) -I$(CPLEX_INCLUDE_PATH) -I$(CONCERT_INCLUDE_PATH) -I$(SCIP_INCLUDE_PATH)
CXXLNFLAGS = -lconcert -lilocplex -lcplex -lm -lpthread -lz -lgmp -lreadline -lscip-3.1.1.linux.x86_64.gnu.opt -llpicpx-3.1.1.linux.x86_64.gnu.opt -lnlpi.cppad-3.1.1.linux.x86_64.gnu.opt -lzimpl.linux.x86_64.gnu.opt

PROG = blk_analysis

all:: $(PROG)

SRC = blocking_analysis.cc \
	taskset_generator.cc \
	main.cc

OBJS = blocking_analysis.o \
	taskset_generator.o \
	main.o

$(OBJS): $(SRC)
	$(CXX) $(CXXFLAGS) -c $(SRC)

$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) $(CXXLN_LIBS_PATH) $(OBJS) -o $(PROG) $(CXXLNFLAGS)

clean::
	rm -rf *.o
	rm blk_analysis
