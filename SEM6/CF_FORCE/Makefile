.PHONY: all run clean

CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall
EXECUTABLE = main
PYTHON_SCRIPT = converterbox.py

all: $(EXECUTABLE)

$(EXECUTABLE): final_optimised.cpp
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) final_optimised.cpp

run: $(EXECUTABLE)
	./$(EXECUTABLE)
	python3 $(PYTHON_SCRIPT)

.PHONY: sim
sim: $(EXECUTABLE)
	./$(EXECUTABLE)

.PHONY: vis
vis:
	python3 $(PYTHON_SCRIPT)