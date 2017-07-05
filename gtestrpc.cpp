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

void
extractCorners(const string& imd_fname,
               vector<ground_coord_type>& gcnrs,
               vector<image_coord_type>&  icnrs)
{
   gcnrs.clear();
   icnrs.clear();
   ifstream istr(imd_fname.c_str());
   long nrows=-1, ncols=-1;
   string key, eq, val;

   do {
      key.clear();
      istr >> key;
      if        (key == "numRows") {
         istr >> eq >> nrows;
      } else if (key == "numColumns") {
         istr >> eq >> ncols;
      }
   } while (key.length() &&
            key != "BEGIN_GROUP");
   getline(istr, key); // clear the rest of BEGIN_GROUP line

   if (nrows < 0 || ncols < 0)
      return; // empty vectors ==> failure

   image_coord_type ip;
   ip.s[0] = 0;
   ip.s[1] = 0;
   icnrs.push_back(ip);
   ip.s[0] = ncols-1;
   icnrs.push_back(ip);
   ip.s[1] = nrows-1;
   icnrs.push_back(ip);
   ip.s[0] = 0;
   icnrs.push_back(ip);
   
   // TBD: why can't I use initializing constructors?
   ground_coord_type gp;// = (cl_double3)(0.0,0.0,0.0);
   gcnrs = vector<ground_coord_type>(4, gp);

   do {
      key.clear();
      val.clear();
      istr >> key >> eq >> val;
      //cout << key << " " << eq << " " << val << endl;
      double v = atof(val.c_str());
      if      (key=="ULLon")
         gcnrs[0].s[0] = v; 
      else if (key=="ULLat")
         gcnrs[0].s[1] = v; 
      else if (key=="ULHAE")
         gcnrs[0].s[2] = v; 
      else if (key=="URLon")
         gcnrs[1].s[0] = v; 
      else if (key=="URLat")
         gcnrs[1].s[1] = v; 
      else if (key=="URHAE")
         gcnrs[1].s[2] = v; 
      else if (key=="LRLon")
         gcnrs[2].s[0] = v; 
      else if (key=="LRLat")
         gcnrs[2].s[1] = v; 
      else if (key=="LRHAE")
         gcnrs[2].s[2] = v; 
      else if (key=="LLLon")
         gcnrs[3].s[0] = v; 
      else if (key=="LLLat")
         gcnrs[3].s[1] = v; 
      else if (key=="LLHAE")
         gcnrs[3].s[2] = v;
   }
   while (key != "END_GROUP" && val.length());

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
      vector <image_coord_type> icnrs;
      extractCorners(base+".IMD", gcnrs, icnrs);
      ASSERT_EQ(4, gcnrs.size());
      ASSERT_EQ(4, icnrs.size());

      ground_coord_type* gc = first(gcnrs);
      image_coord_type ip; //(0,0);
      vector<image_coord_type> g2is(4,ip);
      image_coord_type* ic = first(g2is);
      e = rpc.llh2sl(4, gc, ic);
      EXPECT_SUCCESS(e);

      for (size_t i=0; i<gcnrs.size(); ++i) {
         EXPECT_NEAR(icnrs[i].s[0], g2is[i].s[0], 0.5);
         EXPECT_NEAR(icnrs[i].s[1], g2is[i].s[1], 0.5);
      }
   }
}
