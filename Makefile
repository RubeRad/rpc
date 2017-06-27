CC   = g++
EXE  = rpc
OBJS = rpc.o main.o

$(EXE): $(OBJS)

rpc.o: rpc.cpp rpc.h

main.o: main.cpp rpc.h

clean:
	rm -rf $(EXE) $(OBJS)
