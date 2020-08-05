# compiler

# libraries
ARMA = -l armadillo
LIBS = $(ARMA)

# compiler flags
CXXFLAGS	= -std=c++17 -Ofast -march=native -g -Wall -fopenmp
OBJFILES 	= main.o 
TARGET   	= ../bin/arma_test

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CXX) $(CXXFLAGS) $(LIBS) -o $(TARGET) $(OBJFILES)

clean:
	rm -v $(TARGET) $(OBJFILES)
