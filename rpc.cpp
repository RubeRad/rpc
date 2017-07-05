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
                 ground_coord_type* llh,  // 3n array of lon,lat,hae (deg/m)
                 image_coord_type*  sl)   // preallocated 2n array for samp,line
{
   for (size_t i=0; i<n; ++i) {
      llh2sl_single(&off_scl[0], &coeffs[0],
                    llh[i*3], llh[i*3+1], llh[i*3+2],
                     sl[i*2],  sl[i*2+1]);
   }
   return 0;
}


// TBD turn into an openCL kernel
void
RPC::llh2sl_single(normalizer_type*  off_scl,
                   coefficient_type* coeffs,
                   ground_coord_type lon,
                   ground_coord_type lat,
                   ground_coord_type hae,
                   image_coord_type& samp,
                   image_coord_type& line)
{
   // normalize the input ground point
   image_coord_type x = (image_coord_type)(lon - off_scl[0]) / off_scl[5];
   image_coord_type y = (image_coord_type)(lat - off_scl[1]) / off_scl[6];
   image_coord_type z = (image_coord_type)(hae - off_scl[2]) / off_scl[7];
   image_coord_type xx=x*x, yy=y*y, zz=z*z, xy=x*y, xz=x*z, yz=y*z;
   image_coord_type xxx=xx*x, xxy=xx*y, xxz=xx*z,
                     xyy=x*yy, yyy=yy*y, yyz=yy*z,
                     xzz=x*zz, yzz=y*zz, zzz=zz*z, xyz=x*y*z;
   
   image_coord_type sampn=0, sampd=0, linen=0, lined=0;
   // accumulate in reverse order, adding smallest terms first
   for (int i=19; i>=0; --i) {
      image_coord_type term;
      switch(i) {
      case  0: term =   1; break;
      case  1: term =   x; break;
      case  2: term =   y; break;
      case  3: term =   z; break;
      case  4: term =  xy; break;
      case  5: term =  xz; break;
      case  6: term =  yz; break;
      case  7: term =  xx; break;
      case  8: term =  yy; break;
      case  9: term =  zz; break;
      case 10: term = xyz; break;
      case 11: term = xxx; break;
      case 12: term = xyy; break;
      case 13: term = xzz; break;
      case 14: term = xxy; break;
      case 15: term = yyy; break;
      case 16: term = yzz; break;
      case 17: term = xxz; break;
      case 18: term = yyz; break;
      case 19: term = zzz; break;
      }
      // coeff groups are lineNum, lineDen, sampNum, sampDen as per RPB file
      linen += term * coeffs[i];
      lined += term * coeffs[i+20];
      sampn += term * coeffs[i+40];
      sampd += term * coeffs[i+60];
   }

   image_coord_type soff=(image_coord_type)off_scl[3];
   image_coord_type loff=(image_coord_type)off_scl[4];
   image_coord_type sscl=(image_coord_type)off_scl[8];
   image_coord_type lscl=(image_coord_type)off_scl[9];
   if (sampd == 0 || lined == 0) {
      samp = (image_coord_type)-1e10;
      line = (image_coord_type)-1e10;
   } else {
      samp = sampn / sampd * sscl + soff;
      line = linen / lined * lscl + loff;
   }
}




