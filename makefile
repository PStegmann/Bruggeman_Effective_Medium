CXX=clang++ -std=c++17 -Xpreprocessor -fopenmp
CFLAGS=-c -O0 -g -Wall -Wextra -Wno-unused-parameter -Wc++11-extensions -Ic:/usr/local/Cellar/boost/1.76.0/include/
LIBS=-lm -lomp -lboost_system -lboost_filesystem
SOURCES=get_refractive_index_sand.cpp calcite.cpp ferrum.cpp quartz.cpp silicate.cpp soot.cpp sulphates.cpp water.cpp
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=get_index_exec

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LIBS) -o $@

.c.o:
	$(CXX) $(CFLAGS) $< -o $@
