.PHONY: all run batch sim vis clean

CXX      = g++
CXXFLAGS = -O3 -std=c++17 -Wall
EXECUTABLE = main
SOURCE     = final_optimised.cpp
PYTHON_SCRIPT = convert.py

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $(SOURCE)

# single run — default init (rings) or pass INIT=name
# usage: make run
#        make run INIT=vortex
run: $(EXECUTABLE)
	./$(EXECUTABLE) $(INIT)

# run ALL init conditions sequentially
batch: $(EXECUTABLE)
	python3 convert.py

# run a subset — usage: make subset INITS="rings vortex lattice"
subset: $(EXECUTABLE)
	python3 convert.py $(INITS)

# manually re-convert a folder
# usage: make vis FOLDER=run_rings_N1000_...
vis:
	python3 $(PYTHON_SCRIPT) $(FOLDER)

clean:
	rm -f $(EXECUTABLE)