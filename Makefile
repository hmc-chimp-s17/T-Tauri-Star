CXX = clang++

CXX_FLAGS = -g -std=c++11 -pedantic -Wall -Wextra

all: simulation

# Compile "ttauristar.cpp" to create object file "ttauristar.o"
ttauristar.o: ttauristar.cpp ttauristar.hpp
	$(CXX) -c -o ttauristar.o $(CXX_FLAGS) ttauristar.cpp
	
# Compile "simulation.cpp" to create object file "simulation.o"
simulation.o: simulation.cpp 
	$(CXX) -c -o simulation.o $(CXX_FLAGS) simulation.cpp

# Create executable "simulation" by linking 
simulation: simulation.o ttauristar.o
	$(CXX) -o simulation ttauristar.o simulation.o

clean:
	rm -rf *.o simulation