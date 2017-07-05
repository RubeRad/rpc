#include <gtest/gtest.h>
#include <fstream>

#include "rpc.h"

#define EXPECT_SUCCESS(x) EXPECT_EQ(0,x)
#define ASSERT_SUCCESS(x) ASSERT_EQ(0,x)

using namespace RPC;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;

vector<ground_coord_type>
extractCorners(const string& imd_fname) {
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


TEST(RPC, g2iCorners) {
   vector<string> bases;
   bases.push_back("data/13DEC28032941-M1BS-053950035030_01_P001");
   bases.push_back("data/13DEC28032941-P1BS-053950035030_01_P001");
   bases.push_back("data/13DEC28033039-M1BS-053950035030_01_P001");
   bases.push_back("data/13DEC28033039-P1BS-053950035030_01_P001");

   for (const auto& base : bases) {
      cout << "Testing: " << base << endl;
      RPC::RPC rpc;
      errorType e = rpc.init(base+".RPB");
      ASSERT_SUCCESS(e);

      vector<ground_coord_type> gcnrs;
      gcnrs = extractCorners(base+".IMD");
      ASSERT_SUCCESS(e);

      ground_coord_type* gc = first(gcnrs);
      vector<image_coord_type> icnrs(8,0);
      image_coord_type* ic = first(icnrs);
      e = rpc.llh2sl(4, gc, ic);
      EXPECT_SUCCESS(e);
   }
}
