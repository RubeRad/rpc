#include <gtest/gtest.h>
#include <fstream>
#include <cmath>
#include <random>

#include "rpc.hpp"

#define EXPECT_SUCCESS(x) EXPECT_EQ(0,x)
#define ASSERT_SUCCESS(x) ASSERT_EQ(0,x)

using namespace RPC_NS;
using std::cout;
using std::endl;
using std::vector;
using std::string;
using std::ifstream;

void
extractCorners(const string& imd_fname,
               vector<ground_coord_type>& gcnrs,
               vector<imaged_coord_type>& icnrs)
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

   imaged_coord_type ip;
   ip.x = 0;
   ip.y = 0;
   icnrs.push_back(ip);
   ip.x = ncols-1;
   icnrs.push_back(ip);
   ip.y = nrows-1;
   icnrs.push_back(ip);
   ip.x = 0;
   icnrs.push_back(ip);
   
   ground_coord_type gp;// = (cl_double3)(0.0,0.0,0.0);
   gcnrs = vector<ground_coord_type>(4, gp);

   do {
      key.clear();
      val.clear();
      istr >> key >> eq >> val;
      //cout << key << " " << eq << " " << val << endl;
      double v = atof(val.c_str());
      if      (key=="ULLon")  gcnrs[0].x = v; 
      else if (key=="ULLat")  gcnrs[0].y = v; 
      else if (key=="ULHAE")  gcnrs[0].z = v; 
      else if (key=="URLon")  gcnrs[1].x = v; 
      else if (key=="URLat")  gcnrs[1].y = v; 
      else if (key=="URHAE")  gcnrs[1].z = v; 
      else if (key=="LRLon")  gcnrs[2].x = v; 
      else if (key=="LRLat")  gcnrs[2].y = v; 
      else if (key=="LRHAE")  gcnrs[2].z = v; 
      else if (key=="LLLon")  gcnrs[3].x = v; 
      else if (key=="LLLat")  gcnrs[3].y = v; 
      else if (key=="LLHAE")  gcnrs[3].z = v;
   }
   while (key != "END_GROUP" && val.length());

}

#if 0
cl_double8
g2i_partials_numerical(const RPC::RPC& rpc,
                       const ground_coord_type& gp,
                       double step_scl=1.0)
{
   ground_coord_type gps[4];
   imaged_coord_type ips[4];
   for (int i=0; i<4; ++i)
      gps[i] = gp;
   double dh = 0.00001*step_scl, dv = 1.0*step_scl;
   gps[1].x += dh;
   gps[2].y += dh;
   gps[3].z += dv;
   rpc.llh2sl(4, &gps[0], &ips[0]);

   cl_double8 ansa;
   ansa.x = ips[0].x;
   ansa.y = ips[0].y;
   ansa.s2 = (ips[1].x-ips[0].x)/dh;  ansa.s5 = (ips[1].y-ips[0].y)/dh; // ds/dX, dl/dX
   ansa.s3 = (ips[2].x-ips[0].x)/dh;  ansa.s6 = (ips[2].y-ips[0].y)/dh; // ds/dY, dl/dY
   ansa.s4 = (ips[3].x-ips[0].x)/dv;  ansa.s7 = (ips[3].y-ips[0].y)/dv; // ds/dZ, dl/dZ

   return ansa;
}
#endif

void EXPECT_NEAR_PCT(double shld,
                     double isss,
                     double pct,
                     int index=-1,
                     const string& lbl="")
{
   double tol = fabs(shld) * pct / 100.0;
   EXPECT_NEAR(shld, isss, tol) << "index: " << index << " " << lbl;
}

std::vector<std::string>
test_files(const std::string& ext = "")
{
   std::vector<std::string> bases;
   bases.push_back("data/13DEC28032941-M1BS-053950035030_01_P001");
   bases.push_back("data/13DEC28032941-P1BS-053950035030_01_P001");
   bases.push_back("data/13DEC28033039-M1BS-053950035030_01_P001");
   bases.push_back("data/13DEC28033039-P1BS-053950035030_01_P001");
   if (ext.length())
      for (auto& b : bases)
         b += ext;
   return bases;
}


TEST(RPC, g2iCorners) {
   auto bases = test_files();
   for (const auto& base : bases) {
      cout << "Testing: " << base << endl;
      RPC<double> rpc;
      errorType e = rpc.init(base+".RPB");
      ASSERT_SUCCESS(e);

      vector<ground_coord_type> gcnrs;
      vector<imaged_coord_type> icnrs;
      extractCorners(base+".IMD", gcnrs, icnrs);
      ASSERT_EQ(4, gcnrs.size());
      ASSERT_EQ(4, icnrs.size());

      imaged_coord_type ip;
      vector<imaged_coord_type> g2is;
      for (auto& gp : gcnrs) {
         rpc.g2i(gp, ip);
         g2is.push_back(ip);
      }

      // The ground corners in the IMD file come from the rigorous sensor model
      // The RPC in the RPB file will have some amount of fit error
      // For this data mostly the corners are <0.1, but there is a .10, .11, .13
      const double fit_tol = 0.15;
      for (size_t i=0; i<gcnrs.size(); ++i) {
         EXPECT_NEAR(icnrs[i].x, g2is[i].x, fit_tol) << base << " " << i << " X ";
         EXPECT_NEAR(icnrs[i].y, g2is[i].y, fit_tol) << base << " " << i << " Y ";



#if 0
         cl_double8 sl_part, sl_part_num;
         g2ipartials(&rpc.off_scl[0], &rpc.coeffs[0], gcnrs[i], sl_part);
         EXPECT_NEAR(g2is[i].x, sl_part.x, 1.0e-12);
         EXPECT_NEAR(g2is[i].y, sl_part.y, 1.0e-12);

         sl_part_num = g2i_partials_numerical(rpc, gcnrs[i]);

         double scl = 64;
         for (int it=0; it<10; ++it) {
            cl_double8 dum = g2i_partials_numerical(rpc, gcnrs[i], scl);
            scl /= 2;
            cout << "scl " << scl << "\t" << dum.s[3] << "\t" << sl_part.s[3] << endl;
         }

         EXPECT_NEAR(sl_part_num.s[0], sl_part.s[0], 1.0e-12); // samp
         EXPECT_NEAR(sl_part_num.s[1], sl_part.s[1], 1.0e-12); // line

         EXPECT_NEAR_PCT(sl_part_num.s[2], sl_part.s[2], 1,  2); // ds/dx significant
         EXPECT_NEAR_PCT(sl_part_num.s[3], sl_part.s[3], 10, 3);
         EXPECT_NEAR_PCT(sl_part_num.s[4], sl_part.s[4], 10, 4);
         EXPECT_NEAR_PCT(sl_part_num.s[5], sl_part.s[5], 10, 5);
         EXPECT_NEAR_PCT(sl_part_num.s[6], sl_part.s[6], 1,  6); // dl/dy significant
         EXPECT_NEAR_PCT(sl_part_num.s[7], sl_part.s[7], 10, 7);
         continue;

         double norml = (icnrs[i].x - rpc.off_scl[OFFS]) / rpc.off_scl[SCLS];
         double norms = (icnrs[i].y - rpc.off_scl[OFFL]) / rpc.off_scl[SCLL];
         double normz = (gcnrs[i].z - rpc.off_scl[OFFZ]) / rpc.off_scl[SCLZ];
         double normx, normy;
         i2g_dlt(rpc.coeffs, norml, norms, normz, normx, normy);
         double dltx = normx * rpc.off_scl[SCLX] + rpc.off_scl[OFFX];
         double dlty = normy * rpc.off_scl[SCLY] + rpc.off_scl[OFFY];
         double rtl = ( rpc.coeffs[COEFF_1   ] +
                        rpc.coeffs[COEFF_X   ] * normx +
                        rpc.coeffs[COEFF_Y   ] * normy +
                        rpc.coeffs[COEFF_Z   ] * normz ) /
                      ( rpc.coeffs[COEFF_1+20] +
                        rpc.coeffs[COEFF_X+20] * normx +
                        rpc.coeffs[COEFF_Y+20] * normy +
                        rpc.coeffs[COEFF_Z+20] * normz );
         double rts = ( rpc.coeffs[COEFF_1+40] +
                        rpc.coeffs[COEFF_X+40] * normx +
                        rpc.coeffs[COEFF_Y+40] * normy +
                        rpc.coeffs[COEFF_Z+40] * normz ) /
                      ( rpc.coeffs[COEFF_1+60] +
                        rpc.coeffs[COEFF_X+60] * normx +
                        rpc.coeffs[COEFF_Y+60] * normy +
                        rpc.coeffs[COEFF_Z+60] * normz );
         EXPECT_NEAR(norml, rtl, 1.0e-12);
         EXPECT_NEAR(norms, rts, 1.0e-12);
         EXPECT_NEAR(gcnrs[i].x, dltx, 1.0e-4); // ~ 10m
         EXPECT_NEAR(gcnrs[i].y, dlty, 1.0e-4);
#endif
      } // for i
   } // for base
} // TEST


TEST(RPC, doubleVSfloat) {
   std::mt19937 seed(828067); // ASCII 82 80 67 = R P C
   std::uniform_real_distribution<double> uni(-1.0, 1.0);
   std::vector<ground_coord_type> gps;
   for (int i=0; i<100; ++i)
      gps.emplace_back(uni(seed), uni(seed), uni(seed));

   auto rpbs = test_files(".RPB");
   for (const auto& rpb : rpbs) {
      RPC<double> rpcd;
      RPC<float>  rpcf;
      ASSERT_SUCCESS(rpcd.init(rpb));
      ASSERT_SUCCESS(rpcf.init(rpb));

      for (const auto& gp1 : gps) {
         imagef_coord_type ipf;
         imaged_coord_type ipd;
         ground_coord_type gp = gp1;
         gp.x = gp1.x*rpcd.off_scl[SCLX] + rpcd.off_scl[OFFX];
         gp.y = gp1.y*rpcd.off_scl[SCLY] + rpcd.off_scl[OFFY];
         gp.z = gp1.z*rpcd.off_scl[SCLZ] + rpcd.off_scl[OFFZ];
         
         rpcf.g2i(gp, ipf);
         rpcd.g2i(gp, ipd);
         EXPECT_NEAR(ipd.x, ipf.x, 0.1);
         EXPECT_NEAR(ipd.y, ipf.y, 0.1);

         std::cout << "RPC_D_VS_F," << ipf.x << "," << ipf.y << ","
                   << ipd.x << "," << ipd.y << std::endl;
      }
   }
}


#if 0
void test_partials(RPC& rpc, ground_coord_type gp, const string& lbl) {
   cl_double8 prig, pnum;
   g2ipartials(&rpc.off_scl[0], &rpc.coeffs[0], gp, prig);
   pnum = g2i_partials_numerical(rpc, gp, 0.1);
   EXPECT_NEAR(pnum.s[0], prig.s[0], 1.0e-12); // samp
   EXPECT_NEAR(pnum.s[1], prig.s[1], 1.0e-12); // line

   EXPECT_NEAR_PCT(pnum.s[2], prig.s[2], .1, 2, lbl); // ds/dx significant
   EXPECT_NEAR_PCT(pnum.s[3], prig.s[3],  1, 3, lbl);
   EXPECT_NEAR_PCT(pnum.s[4], prig.s[4],  1, 4, lbl);
   EXPECT_NEAR_PCT(pnum.s[5], prig.s[5],  1, 5, lbl);
   EXPECT_NEAR_PCT(pnum.s[6], prig.s[6], .1, 6, lbl); // dl/dy significant
   EXPECT_NEAR_PCT(pnum.s[7], prig.s[7],  1, 7, lbl);
}
   

TEST(RPC, partials) {
   string base = "data/13DEC28032941-M1BS-053950035030_01_P001";
   RPC::RPC rpc;
   errorType e = rpc.init(base+".RPB");
   ASSERT_SUCCESS(e);

   ground_coord_type gp0, gpx, gpy, gpz, gpxy, gpxz, gpyz, gpxyz;
   gp0.x = rpc.off_scl[0];
   gp0.y = rpc.off_scl[1];
   gp0.z = rpc.off_scl[2];
   test_partials(rpc, gp0, "0");

   gpx = gp0;
   gpx.x += 0.001;
   test_partials(rpc, gpx, "X");

   gpy = gp0;
   gpy.y += 0.001;
   test_partials(rpc, gpy, "Y");
   
   gpz = gp0;
   gpz.z += 100;
   test_partials(rpc, gpz, "Z");

   gpxy = gp0;
   gpxy.x += 0.001;
   gpxy.y += 0.001;
   test_partials(rpc, gpxy, "XY");   

   gpxz= gp0;
   gpxz.x += 0.001;
   gpxz.z += 100;
   test_partials(rpc, gpxz, "XZ");   

   gpyz= gp0;
   gpyz.y += 0.001;
   gpyz.z += 100;
   test_partials(rpc, gpyz, "YZ");   

   gpxyz= gp0;
   gpxyz.x += 0.001;
   gpxyz.y += 0.001;
   gpxyz.z += 100;
   test_partials(rpc, gpxyz, "XYZ");   
}

#endif
