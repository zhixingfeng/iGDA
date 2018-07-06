CXX = g++
CC = gcc
INCLUDES = -I include -I tools/boost/include -I tools -I tools/stxxl/build/local/stxxl/include -I tools/seqan/seqan/include
LIBS = -L tools/stxxl/build/local/stxxl/lib
CXX_FLAGS = -pthread -std=c++14 -w -O3 -lstxxl -lz

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
	$(wildcard tools/ssw/*.cpp) \
    $(wildcard tools/prob/*.cpp)

SRC_C = $(wildcard tools/ssw/*.c)

OBJ_CXX = $(SRC_CXX:.cpp=.o)

OBJ_C = $(SRC_C:.c=.o)

all: mkbin igda rmobj
	
.PHONY: mkbin
mkbin:
	mkdir -p bin
	
igda: $(OBJ_CXX) $(OBJ_C) 
	$(CXX) -o bin/igda $^ $(LIBS) $(CXX_FLAGS)

%.o: %.cpp
	$(CXX) $(INCLUDES) -c $< -o $@ $(CXX_FLAGS)

%.o: %.c
	$(CC) $(INCLUDES) -c $< -o $@ 
rmobj:
	rm -f $(OBJ_CXX)

.PHONY: clean
clean:
	rm -rf bin

