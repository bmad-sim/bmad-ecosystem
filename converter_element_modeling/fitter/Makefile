CXX=g++
CXXFLAGS=-Wall -Wextra -Wpedantic -Weffc++ -std=c++17
LD=g++
LDFLAGS=-lgsl -lgslcblas -lm

fitter: main.o read_data.o cauchy.o meta_fit.o
	$(LD) $^ -o $@ $(LDFLAGS)

newmain.o: main.cpp read_data.hpp cauchy.hpp

read_data.o: read_data.cpp read_data.hpp cauchy.hpp

cauchy.o: cauchy.cpp cauchy.hpp

meta_fit.o: meta_fit.cpp meta_fit.hpp

.PHONY: clean
clean:
	-rm *.o
	-rm fitter
