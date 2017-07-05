CPPFLAGS = -isystem $(GTEST_DIR)/include
CXXFLAGS = -g -std=c++11 -Wall -Wextra -pthread
OBJS     = rpc.o gtestrpc.o

all: gtestrpc

rpc.o: rpc.cpp rpc.h

gtestrpc.o: gtestrpc.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(GTEST_DIR)/include -c $^

gtestrpc: rpc.o gtestrpc.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -lpthread $^ -o $@ $(GTEST_DIR)/make/gtest_main.a



clean:
	rm -rf gtestrpc $(OBJS)
