PREFIX = .
CXX = g++
INCLUDES = -I include -I tools/boost/include -I tools -I tools/seqan/seqan/include -I tools/stxxl/build/local/stxxl/include
LIBS = -L tools/stxxl/build/local/stxxl/lib
CXX_FLAGS = -pthread -std=c++14 -w -O2 -lz -lstxxl

SRC_CXX = $(wildcard src/*.cpp)\
	$(wildcard src/modules/aligncoder/*.cpp)\
	$(wildcard src/modules/alignreader/*.cpp)\
	$(wildcard src/modules/dforest/*.cpp)\
	$(wildcard src/modules/errormodel/*.cpp)\
	$(wildcard src/modules/hclust/*.cpp)\
	$(wildcard src/modules/sclust/*.cpp)\
	$(wildcard src/modules/assemble/*.cpp)\
    $(wildcard src/modules/detectsingle/*.cpp)\
    $(wildcard src/modules/rsm/*.cpp)\
	$(wildcard test/*.cpp)\
    $(wildcard tools/prob/*.cpp)\
	$(wildcard src/misc/*.cpp)


OBJ_CXX = $(SRC_CXX:.cpp=.o)

all: mkbin igda rmobj
	
.PHONY: mkbin
mkbin:
	mkdir -p $(PREFIX)/bin
	
igda: $(OBJ_CXX)
	$(CXX) -o $(PREFIX)/bin/igda $^ $(LIBS) $(CXX_FLAGS)

%.o: %.cpp
	$(CXX) $(INCLUDES) -c $< -o $@ $(CXX_FLAGS)

rmobj:
	rm -f $(OBJ_CXX)

.PHONY: clean
clean:
	rm -rf $(PREFIX)/bin

