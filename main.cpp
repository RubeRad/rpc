#include "rpc.h"
#include <iostream>

using namespace RPC;
using std::cout;
using std::endl;
using std::vector;

int
main(int argc, char** argv) {
   RPC::RPC rpc;
   errorType e = -1;
   if (argc >= 2)
      e = rpc.init(argv[1]);
   if (IS_OK(e))
      cout << "RPB file OK\n";
   else
      cout << "RPB file not OK\n";

   vector<ground_coord_type> cnrs;
   if (argc >= 3)
      cnrs = extractCorners(argv[2]);
   if (cnrs.empty())
      cout << "IMD file not OK\n";
   else 
      cout << "IMD file OK\n";

   double* c = first(cnrs);
   for (size_t i=0; i<cnrs.size(); i++) {
      cout << i << "\t" << c[i] << endl;
   }
}
