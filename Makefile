CXX = g++

PROG = taskgen

all:: $(PROG) 

SRC = BlockingAnalysis.cpp \
	TaskSetGenerator.cpp \
	TestTaskSetGen.cpp

OBJS = BlockingAnalysis.o \
	TaskSetGenerator.o \
	TestTaskSetGen.o

$(OBJS): $(SRC)
	$(CXX) -c $(SRC)

$(PROG): $(OBJS)
	$(CXX) $(OBJS) -o $(PROG)

clean::
	rm -rf *.o
	rm taskgen
