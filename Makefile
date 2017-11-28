CXX = g++
CC = gcc
INCLUDES = -I include -I tools/boost/include -I tools
CXX_FLAGS = -pthread -std=c++11 -w -O3 

SRC_CXX = $(wildcard src/*.cpp) \
	$(wildcard src/modules/aligncoder/*.cpp)\
	$(wildcard src/modules/alignreader/*.cpp)\
	$(wildcard src/modules/dforest/*.cpp)\
    	$(wildcard src/modules/errormodel/*.cpp)\
	$(wildcard src/modules/hclust/*.cpp)\
	$(wildcard src/modules/sclust/*.cpp)\
	$(wildcard src/modules/assemble/*.cpp)\
	$(wildcard test/*.cpp) \
	$(wildcard tools/boost/src/filesystem/*.cpp) \
	$(wildcard tools/boost/src/system/*.cpp) \
	$(wildcard tools/ssw/*.cpp) 

SRC_C = $(wildcard tools/ssw/*.c)

OBJ_CXX = $(SRC_CXX:.cpp=.o)

OBJ_C = $(SRC_C:.c=.o)

all: mkbin igda rmobj
	
.PHONY: mkbin
mkbin:
	mkdir -p bin
	
igda: $(OBJ_CXX) $(OBJ_C) 
	$(CXX) -o bin/igda $^ $(CXX_FLAGS)

%.o: %.cpp
	$(CXX) $(INCLUDES) -c $< -o $@ $(CXX_FLAGS)

%.o: %.c
	$(CC) $(INCLUDES) -c $< -o $@ 
rmobj:
	rm -f $(OBJ_CXX)

.PHONY: clean
clean:
	rm -rf bin

