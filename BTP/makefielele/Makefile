.PHONY: all run clean sim vis

CXX      = g++
CXXFLAGS = -O3 -std=c++17 -Wall
EXECUTABLE = main
SOURCE     = new_final.cpp
PYTHON_SCRIPT = convert.py

all: $(EXECUTABLE)

$(EXECUTABLE): $(SOURCE)
	$(CXX) $(CXXFLAGS) -o $(EXECUTABLE) $(SOURCE)

# runs sim — conversion happens automatically inside the exe
run: $(EXECUTABLE)
	./$(EXECUTABLE)

# same as run, explicit alias
sim: $(EXECUTABLE)
	./$(EXECUTABLE)

# manually re-run conversion on a specific folder
# usage: make vis FOLDER=run_N500_rho0.5_b10_a10_Dr0.1_dt0.001
vis:
	python3 $(PYTHON_SCRIPT) $(FOLDER)

clean:
	rm -f $(EXECUTABLE)