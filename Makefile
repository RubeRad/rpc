#LD       = g++
CC       = g++ # why doesn't LD=g++ work?
CPPFLAGS = -g
EXE      = rpc
OBJS     = rpc.o main.o

all: $(EXE)

$(EXE): $(OBJS)

rpc.o: rpc.cpp rpc.h

main.o: main.cpp rpc.h

clean:
	rm -rf $(EXE) $(OBJS)
