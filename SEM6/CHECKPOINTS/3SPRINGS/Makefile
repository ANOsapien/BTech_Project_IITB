.PHONY: all run clean

CXX = g++
CXXFLAGS = -O3 -std=c++17 -Wall
EXECUTABLE = main
PYTHON_SCRIPT = converter_2.py

all: $(EXECUTABLE)

$(EXECUTABLE): VE_CF_SPRING.cpp
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) VE_CF_SPRING.cpp

run: $(EXECUTABLE)
	./$(EXECUTABLE)
	python3 $(PYTHON_SCRIPT)

.PHONY: sim
sim: $(EXECUTABLE)
	./$(EXECUTABLE)

.PHONY: vis
vis:
	python3 $(PYTHON_SCRIPT)