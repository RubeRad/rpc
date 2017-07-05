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

errorType
RPC::RPC::llh2sl(size_t n,
                 ground_coord_type* llh, // n-array of lon,lat,hae (deg/m) triplets
                 image_coord_type*  sl)  // preallocated n-array for samp,line duplets
{
   for (size_t i=0; i<n; ++i) {
      // TBD: Why can't I use cl_double3.x, .y, .z?
      llh2sl_single(&off_scl[0], &coeffs[0],
                    llh[i].s[0], llh[i].s[1], llh[i].s[2],
                     sl[i].s[0],  sl[i].s[1]);
   }
   return 0;
}


// TBD turn into an openCL kernel
void
RPC::llh2sl_single(normalizer_type*  off_scl,
                   coefficient_type* coeffs,
                   double  lon,
                   double  lat,
                   double  hae,
                   double& samp,
                   double& line)
{
   // normalize the input ground point
   double x = (lon - off_scl[0]) / off_scl[5];
   double y = (lat - off_scl[1]) / off_scl[6];
   double z = (hae - off_scl[2]) / off_scl[7];
   double xx=x*x, yy=y*y, zz=z*z, xy=x*y, xz=x*z, yz=y*z;
   
   double sampn=0, sampd=0, linen=0, lined=0;
   // accumulate in reverse order, adding smallest terms first
   for (int i=19; i>=0; --i) {
      double term;
      switch(i) {
      case  0: term =    1; break;
      case  1: term =    x; break;
      case  2: term =    y; break;
      case  3: term =    z; break;
      case  4: term =   xy; break;
      case  5: term =   xz; break;
      case  6: term =   yz; break;
      case  7: term =   xx; break;
      case  8: term =   yy; break;
      case  9: term =   zz; break;
      case 10: term = xy*z; break;
      case 11: term = xx*x; break;
      case 12: term = xy*y; break;
      case 13: term = xz*z; break;
      case 14: term = xx*y; break;
      case 15: term = yy*y; break;
      case 16: term = yz*z; break;
      case 17: term = xx*z; break;
      case 18: term = yy*z; break;
      case 19: term = zz*z; break;
      }
      // coeff groups are lineNum, lineDen, sampNum, sampDen as per RPB file
      linen += term * coeffs[i];
      lined += term * coeffs[i+20];
      sampn += term * coeffs[i+40];
      sampd += term * coeffs[i+60];
   }

   double soff = off_scl[3];
   double loff = off_scl[4];
   double sscl = off_scl[8];
   double lscl = off_scl[9];
   if (sampd == 0 || lined == 0) {
      samp = -1e10;
      line = -1e10;
   } else {
      samp = sampn / sampd * sscl + soff;
      line = linen / lined * lscl + loff;
   }
}




