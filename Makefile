CXX = g++
OPTFLAGS = -O3
BASEFLAGS = -g -std=c++2a -Wall -I./include -D_GLIBCXX_PARALLEL
DEFAULTFLAGS = $(BASEFLAGS) -fopenmp
MACFLAGS = $(BASEFLAGS) -Xpreprocessor -fopenmp -lomp $(OPTFLAGS)
CXXFLAGS = $(DEFAULTFLAGS) $(OPTFLAGS)

vpath %.h include
vpath %.cpp src

simulatePaging: simulation.o LruSizesSim.o OSTree.o IncrementAndKill.o
	$(CXX) $(CXXFLAGS) $^ -o simulatePaging

simulation.o: include/LruSizesSim.h include/params.h include/OSTree.h include/CacheSim.h
OSTree.o: include/OSTree.h

%.o: %.cpp %.h
	$(CXX) $(CXXFLAGS) -c $< -o $@

.PHONY: clean
clean:
	rm -f *.o simulatePaging

.PHONY: debug
debug: CXXFLAGS = $(DEFAULTFLAGS)
debug: simulatePaging

.PHONY: mac
mac: CXXFLAGS = $(MACFLAGS)
mac: simulatePaging