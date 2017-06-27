#include "rpc.h"
#include <iostream>

using namespace RPC;

int
main(int argc, char** argv) {
   RPC::RPC rpc;
   errorType e = -1;
   if (argc >= 2)
      e = rpc.init(argv[1]);
   if (IS_OK(e))
      std::cout << "OK\n";
   else
      std::cout << "Not OK\n";
}
