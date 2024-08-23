CXX      = g++
CXXFLAGS = -std=c++17 -O3
PROGRAM  = gdp

all : program

program : main.o algorithm.o util.o test.o structures.o
	$(CXX) $^ -o $(PROGRAM)

main.o : main.cpp util.hpp test.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test.o : test.cpp test.hpp util.hpp algorithm.hpp structures.hpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

algorithm.o : algorithm.cpp algorithm.hpp structures.hpp 
	$(CXX) $(CXXFLAGS) -c $< -o $@

structures.o : structures.cpp structures.hpp util.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

util.o : util.cpp util.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@


.PHONY: clean
clean:
	rm *.o $(PROGRAM)