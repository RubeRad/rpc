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

   vector<ground_coord_type> gcnrs;
   if (argc >= 3)
      gcnrs = extractCorners(argv[2]);
   if (gcnrs.empty())
      cout << "IMD file not OK\n";
   else 
      cout << "IMD file OK\n";

   ground_coord_type* gc = first(gcnrs);
   vector<image_coord_type> icnrs(8,0);
   image_coord_type* ic = first(icnrs);
   e = rpc.llh2sl(4, gc, ic);
   for (size_t i=0; i<icnrs.size(); ++i) {
      cout << icnrs[i] << endl;
   }
   
}
