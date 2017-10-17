CXX = clang++

CXX_FLAGS = -g -std=c++11 -pedantic -Wall -Wextra

all: simulation

# Compile "simulation.cpp" to create object file "simulation.o"
simulation.o: simulation.cpp
	$(CXX) -c -o simulation.o $(CXX_FLAGS) simulation.cpp

# Create executable "simulation" by linking 
simulation: simulation.o
	$(CXX) -o simulation simulation.o