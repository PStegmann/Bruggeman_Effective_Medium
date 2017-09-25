CXX=clang++
CFLAGS=-c -O2 -Wall -Wextra -Wno-unused-parameter -Wc++11-extensions -Ic:/usr/local/Cellar/boost/1.63.0/include/
LIBS=-lm -lboost_system -lboost_filesystem
SOURCES=get_refractive_index_sand.cpp calcite.cpp ferrum.cpp quartz.cpp silicate.cpp soot.cpp sulphates.cpp water.cpp
OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=get_index_exec

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LIBS) -o $@

.c.o:
	$(CXX) $(CFLAGS) $< -o $@
