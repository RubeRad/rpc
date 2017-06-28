#include "rpc.h"

#include <stdlib.h>
#include <iostream>
#include <fstream>



using namespace RPC;
using std::string;
using std::istream;
using std::ifstream;
using std::cout;
using std::endl;
using std::vector;


errorType
RPC::RPC::init(istream& istr)
{
   string key, eq, val;
   for (int i=0; i<6; ++i) { // skip past errRand
      istr >> key >> eq >> val;
      //cout << key << " " << eq << " " << val << endl;
   }

   for (int i=0; i<10; ++i) {
      istr >> key >> eq >> val;
      double v = atof(val.c_str());
      //cout << key << " " << eq << " " << val << ": " << v << endl;
      if      (key ==   "longOffset")  off_scl[0] = v;
      else if (key ==    "latOffset")  off_scl[1] = v;
      else if (key == "heightOffset")  off_scl[2] = v;
      else if (key ==   "sampOffset")  off_scl[3] = v;
      else if (key ==   "lineOffset")  off_scl[4] = v;
      else if (key ==   "longScale")   off_scl[5] = v;
      else if (key ==    "latScale")   off_scl[6] = v;
      else if (key == "heightScale")   off_scl[7] = v;
      else if (key ==   "sampScale")   off_scl[8] = v;
      else if (key ==   "lineScale")   off_scl[9] = v;
      else {
         // unrecognized key!
         return -1;
      }
   }

   int ci=0;
   for (int j=0; j<4; ++j) {
      istr >> key >> eq >> val; // lineNumCoef = (
      //cout << key << " " << eq << " " << val << endl;
      for (int i=0; i<20; ++i) {
         istr >> val;
         //cout << val << endl;
         coeffs[ci++] = atof(val.c_str());
      }
   }

   istr >> key;
   //cout << key << endl;
   if (key != "END_GROUP") {
      // something went wrong!
      return -1;
   }

   return 0;
}

errorType
RPC::RPC::init(const string& fname)
{
   //cout << "fname is " << fname << endl;
   ifstream istr(fname.c_str());
   return init(istr);
}



vector<ground_coord_type>
RPC::extractCorners(const string& imd_fname) {
   vector<ground_coord_type> cnrs(12, 0.0);
   ifstream istr(imd_fname.c_str());

   string line, key, eq, val;
   do {
      line.clear();
      getline(istr, line);
      //cout << "line: " << line << endl;
   } while (line.length() &&
            line.find("BEGIN_GROUP") == std::string::npos);
   
   do {
      key.clear();
      val.clear();
      istr >> key >> eq >> val;
      //cout << key << " " << eq << " " << val << endl;
      double v = atof(val.c_str());
      if      (key=="ULLon")  cnrs[0]  = v; 
      else if (key=="ULLat")  cnrs[1]  = v; 
      else if (key=="ULHAE")  cnrs[2]  = v; 
      else if (key=="URLon")  cnrs[3]  = v; 
      else if (key=="URLat")  cnrs[4]  = v; 
      else if (key=="URHAE")  cnrs[5]  = v; 
      else if (key=="LRLon")  cnrs[6]  = v; 
      else if (key=="LRLat")  cnrs[7]  = v; 
      else if (key=="LRHAE")  cnrs[8]  = v; 
      else if (key=="LLLon")  cnrs[9]  = v; 
      else if (key=="LLLat")  cnrs[10] = v; 
      else if (key=="LLHAE")  cnrs[11] = v;
   }
   while (key != "END_GROUP" && val.length());

   return cnrs;
}
