CXX = g++ -std=c++14
DBG = -g -w
OPT = -Ofast -DNDEBUG
VALGRIND = -g -DNDEBUG

TARGET_A = RXL

OPTIONS = -lnetworkit -lboost_serialization -lboost_program_options -lboost_system -lboost_filesystem -fopenmp -lboost_timer -lpython3.8
INCLUDEPATH = /usr/local/include/ -I/usr/include/valgrind -I/usr/include/python3.8 -I//usr/local/lib/python3.8/dist-packages/numpy/core/include/
PATHLIB = /usr/local/include/
SOURCES_A = $(TARGET_A).cpp Auxiliary.cpp  Labeling.cpp Labeling_Tools.cpp UpdateData.cpp LabelEntry.cpp Tree.cpp SamPG.cpp Dijkstra.cpp BestRandom.cpp InputOutput.cpp


debug:
	$(CXX) $(DBG) -o $(TARGET_A) -I$(INCLUDEPATH) -L$(PATHLIB) $(SOURCES_A) $(OPTIONS)

all:
	$(CXX) $(DBG) -o $(TARGET_A) -I$(INCLUDEPATH) -L$(PATHLIB) $(SOURCES_A) $(OPTIONS)

release:
	$(CXX) $(OPT) -o $(TARGET_A) -I$(INCLUDEPATH) -L$(PATHLIB) $(SOURCES_A) $(OPTIONS)

valgrind:
	$(CXX) $(VALGRIND) -o $(TARGET_A) -I$(INCLUDEPATH) -L$(PATHLIB) $(SOURCES_A) $(OPTIONS)

clean:
	rm -rf $(TARGET_A)
