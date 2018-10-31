CC=g++
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=main.cpp database.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=BAMBI2b

all: $(SOURCES) $(EXECUTABLE)
	
$(EXECUTABLE): $(OBJECTS) 
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@