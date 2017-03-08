CC = g++
INCLUDES = -I include -I tools/boost/include -I tools
CXX_FLAGS = -pthread -std=c++11 -w -O3 

SRC = $(wildcard src/*.cpp) \
	$(wildcard src/modules/aligncoder/*.cpp)\
	$(wildcard src/modules/alignreader/*.cpp)\
	$(wildcard src/modules/dforest/*.cpp)\
	$(wildcard test/*.cpp) \
	$(wildcard tools/boost/src/filesystem/*.cpp) \
	$(wildcard tools/boost/src/system/*.cpp) 

OBJ = $(SRC:.cpp=.o)

all: mkbin igda rmobj
	
.PHONY: mkbin
mkbin:
	mkdir -p bin
	
igda: $(OBJ) 
	$(CC) -o bin/igda $^ $(CXX_FLAGS)

%.o: %.cpp
	$(CC) $(INCLUDES) -c $< -o $@ $(CXX_FLAGS)

rmobj:
	rm -f $(OBJ)

.PHONY: clean
clean:
	rm -rf bin

