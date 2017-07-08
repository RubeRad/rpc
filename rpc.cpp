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
      llh2sl_single(&off_scl[0], &coeffs[0],
                    llh[i].x, llh[i].y, llh[i].z,
                     sl[i].x,  sl[i].y);
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
      case COEFF_1:   term =    1; break;
      case COEFF_X:   term =    x; break;
      case COEFF_Y:   term =    y; break;
      case COEFF_Z:   term =    z; break;
      case COEFF_XY:  term =   xy; break;
      case COEFF_XZ:  term =   xz; break;
      case COEFF_YZ:  term =   yz; break;
      case COEFF_XX:  term =   xx; break;
      case COEFF_YY:  term =   yy; break;
      case COEFF_ZZ:  term =   zz; break;
      case COEFF_XYZ: term = xy*z; break;
      case COEFF_XXX: term = xx*x; break;
      case COEFF_XYY: term = xy*y; break;
      case COEFF_XZZ: term = xz*z; break;
      case COEFF_XXY: term = xx*y; break;
      case COEFF_YYY: term = yy*y; break;
      case COEFF_YZZ: term = yz*z; break;
      case COEFF_XXZ: term = xx*z; break;
      case COEFF_YYZ: term = yy*z; break;
      case COEFF_ZZZ: term = zz*z; break;
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




