CPPFLAGS = -isystem $(GTEST_DIR)/include
#CXXFLAGS = -g -std=c++11 -Wall -Wextra -pthread -U__STRICT_ANSI__
CXXFLAGS = -g -Wall -Wextra -pthread -U__STRICT_ANSI__ -fsanitize=address
OBJS     = rpc.o gtestrpc.o

all: gtestrpc

gtestrpc.o: gtestrpc.cpp rpc.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(GTEST_DIR)/include -c $^

gtestrpc: gtestrpc.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -lgtest_main -lgtest -lpthread $^ -o $@ 



clean:
	rm -rf gtestrpc $(OBJS)
