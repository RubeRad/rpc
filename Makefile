CPPFLAGS = -isystem $(GTEST_DIR)/include
#CXXFLAGS = -g -std=c++11 -Wall -Wextra -pthread -U__STRICT_ANSI__
CXXFLAGS = -g -Wall -Wextra -pthread -U__STRICT_ANSI__ -fsanitize=address
LDFLAGS = -lgtest_main -lgtest -lpthread -lopencv_core -lopencv_calib3d
OBJS     = rpc.o gtestrpc.o

all: gtestrpc

rpc_ocv.o: rpc_ocv.cpp rpc.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I/usr/include/opencv4 -c $^

gtestrpc.o: gtestrpc.cpp rpc.hpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -I$(GTEST_DIR)/include -c $^

gtestrpc: rpc_ocv.o gtestrpc.o
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) $^ -o $@





clean:
	rm -rf gtestrpc $(OBJS)
