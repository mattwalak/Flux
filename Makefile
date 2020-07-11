all: perlin

SOURCES    = perlin.cpp
OBJECTS    = $(SOURCES:.cpp=.o)

.cpp.o:
	g++ -w -O3 -c $< -o $@

perlin: perlin.o
	g++ $(OBJECTS) $(LDFLAGS) -o $@
	rm perlin.o

clean:
	rm -f *.o

