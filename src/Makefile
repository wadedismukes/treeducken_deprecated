# compiler
CXX = g++

# control variables
CXXFLAGS = -Wall -g
VARNAME = value

objects = Treeducken.o SpeciesTree.o Simulator.o GeneTree.o LocusTree.o MbRandom.o Tree.o Engine.o

install: $(objects) 
	$(CXX) -o ../treeducken $(objects)

Treeducken.o: Treeducken.cpp SpeciesTree.h Simulator.h GeneTree.h LocusTree.h MbRandom.h Tree.h Engine.h
	$(CXX) $(CXXFLAGS) -c Treeducken.cpp

SpeciesTree.o: SpeciesTree.h Tree.h
	$(CXX) $(CXXFLAGS) -c SpeciesTree.cpp

Simulator.o: Simulator.h GeneTree.h
	$(CXX) $(CXXFLAGS) -c Simulator.cpp

GeneTree.o: GeneTree.h LocusTree.h
	$(CXX) $(CXXFLAGS) -c GeneTree.cpp

LocusTree.o: LocusTree.h 
	$(CXX) $(CXXFLAGS) -c LocusTree.cpp

MbRandom.o: MbRandom.h
	$(CXX) $(CXXFLAGS) -c MbRandom.cpp

Tree.o: Tree.h MbRandom.h
	$(CXX) $(CXXFLAGS) -c Tree.cpp

Engine.o: Engine.h Simulator.h
	$(CXX) $(CXXFLAGS) -c Engine.cpp
 
.PHONY : clean 
clean: 
	-rm ../treeducken $(objects)
