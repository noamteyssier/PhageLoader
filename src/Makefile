# compiler

# libraries
ARMA = -l armadillo
BOOST_PO = -l boost_program_options
LIBS = $(ARMA) $(BOOST_PO)

# compiler flags
CXXFLAGS	= -std=c++17 -Ofast -march=native -g -Wall -fopenmp

PREP_OBJFILES 	= preprocess.o
PREP_TARGET   	= ../bin/preprocess

ENRICHMENT_OBJFILES 	= enrichment.o
ENRICHMENT_TARGET   	= ../bin/enrichment

all: $(PREP_TARGET) $(ENRICHMENT_TARGET)

$(PREP_TARGET): $(PREP_OBJFILES)
	$(CXX) $(CXXFLAGS) $(LIBS) -o $(PREP_TARGET) $(PREP_OBJFILES)

$(ENRICHMENT_TARGET): $(ENRICHMENT_OBJFILES)
	$(CXX) $(CXXFLAGS) $(LIBS) -o $(ENRICHMENT_TARGET) $(ENRICHMENT_OBJFILES)

clean:
	rm -v $(PREP_TARGET) $(PREP_OBJFILES)
	rm -v $(ENRICHMENT_TARGET) $(ENRICHMENT_OBJFILES)
