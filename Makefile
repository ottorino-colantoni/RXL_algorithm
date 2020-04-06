CXX = g++ -std=c++14
DBG = -g -w
OPT = -Ofast -DNDEBUG
VALGRIND = -g -DNDEBUG

TARGET_A = RXL

OPTIONS = -lnetworkit -lboost_serialization -lboost_program_options -lboost_system -lboost_filesystem -fopenmp -lboost_timer
INCLUDEPATH = /home/f4b3r/networkit/include/ -I/usr/include/valgrind
PATHLIB = /home/f4b3r/networkit/include
SOURCES_A = $(TARGET_A).cpp Auxiliary.cpp  Labeling.cpp Labeling_Tools.cpp UpdateData.cpp LabelEntry.cpp Tree.cpp SamPG.cpp Dijkstra.cpp

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
