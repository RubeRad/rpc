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
      if      (key ==   "longOffset")  off_scl[OFFX] = v;
      else if (key ==    "latOffset")  off_scl[OFFY] = v;
      else if (key == "heightOffset")  off_scl[OFFZ] = v;
      else if (key ==   "sampOffset")  off_scl[OFFS] = v;
      else if (key ==   "lineOffset")  off_scl[OFFL] = v;
      else if (key ==   "longScale")   off_scl[SCLX] = v;
      else if (key ==    "latScale")   off_scl[SCLY] = v;
      else if (key == "heightScale")   off_scl[SCLZ] = v;
      else if (key ==   "sampScale")   off_scl[SCLS] = v;
      else if (key ==   "lineScale")   off_scl[SCLL] = v;
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
   double x = (lon - off_scl[OFFX]) / off_scl[SCLX];
   double y = (lat - off_scl[OFFY]) / off_scl[SCLY];
   double z = (hae - off_scl[OFFZ]) / off_scl[SCLZ];
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

   if (sampd == 0 || lined == 0) {
      samp = -1e10;
      line = -1e10;
   } else {
      samp = sampn / sampd * off_scl[SCLS] + off_scl[OFFS];
      line = linen / lined * off_scl[SCLL] + off_scl[OFFL];
   }
}



void
RPC::g2ipartials(normalizer_type* off_scl,
                 coefficient_type* coeffs,
                 const cl_double3& llh,
                 cl_double8& sl_part)
{
   // normalize the input ground point
   double x = (llh.x - off_scl[OFFX]) / off_scl[SCLX];
   double y = (llh.y - off_scl[OFFY]) / off_scl[SCLY];
   double z = (llh.z - off_scl[OFFZ]) / off_scl[SCLZ];
   double xx=x*x, yy=y*y, zz=z*z, xy=x*y, xz=x*z, yz=y*z;
   double ln=0, ld=0, sn=0, sd=0;   // line/samp numer/denominator accumulators
   // accumulators for line and sample num/den xyz partials
   double lndx=0, lndy=0, lndz=0, lddx=0, lddy=0, lddz=0;
   double sndx=0, sndy=0, sndz=0, sddx=0, sddy=0, sddz=0;

   // accumulate in reverse order, adding smallest terms first
   for (int i=19; i>=0; --i) {
     double t, tdx, tdy, tdz; // terms x^?y^?z^?, and their x,y,z partials
     switch(i) {
       case COEFF_1:   t=1;    tdx=0;    tdy=0;    tdz=0;    break;
       case COEFF_X:   t=x;    tdx=1;    tdy=0;    tdz=1;    break;
       case COEFF_Y:   t=y;    tdx=0;    tdy=1;    tdz=1;    break;
       case COEFF_Z:   t=z;    tdx=0;    tdy=0;    tdz=1;    break;
       case COEFF_XY:  t=xy;   tdx=y;    tdy=x;    tdz=0;    break;
       case COEFF_XZ:  t=xz;   tdx=z;    tdy=0;    tdz=x;    break;
       case COEFF_YZ:  t=yz;   tdx=0;    tdy=z;    tdz=y;    break;
       case COEFF_XX:  t=xx;   tdx=2*x;  tdy=0;    tdz=0;    break;
       case COEFF_YY:  t=yy;   tdx=0;    tdy=2*y;  tdz=0;    break;
       case COEFF_ZZ:  t=zz;   tdx=0;    tdy=0;    tdz=2*z;  break;
       case COEFF_XYZ: t=xy*z; tdx=yz;   tdy=xz;   tdz=xy;   break;
       case COEFF_XXX: t=xx*x; tdx=3*xx; tdy=0;    tdz=0;    break;
       case COEFF_XYY: t=xy*y; tdx=yy;   tdy=2*xy; tdz=0;    break;
       case COEFF_XZZ: t=xz*z; tdx=zz;   tdy=0;    tdz=2*xz; break;
       case COEFF_XXY: t=xx*y; tdx=2*xy; tdy=xx;   tdz=0;    break;
       case COEFF_YYY: t=yy*y; tdx=0;    tdy=3*yy; tdz=0;    break;
       case COEFF_YZZ: t=yz*z; tdx=0;    tdy=zz;   tdz=2*yz; break;
       case COEFF_XXZ: t=xx*z; tdx=2*xz; tdy=0;    tdz=xx;   break;
       case COEFF_YYZ: t=yy*z; tdx=0;    tdy=2*yz; tdz=yy;   break;
       case COEFF_ZZZ: t=zz*z; tdx=0;    tdy=0;    tdz=3*zz; break;
     }
     // coeff groups are lineNum, lineDen, sampNum, sampDen as per RPB file
     double cln=coeffs[i], cld=coeffs[i+20], csn=coeffs[i+40], csd=coeffs[i+60];
     ln+=t*cln; lndx+=tdx*cln; lndy+=tdy*cln; lndz+=tdz*cln;
     ld+=t*cld; lddx+=tdx*cld; lddy+=tdy*cld; lddz+=tdz*cld;
     sn+=t*csn; sndx+=tdx*csn; sndy+=tdy*csn; sndz+=tdz*csn;
     sd+=t*csd; sddx+=tdx*csd; sddy+=tdy*csd; sddz+=tdz*csd;
   }
     
   if (ld == 0 || sd == 0) {
      for (int i=0; i<8; ++i)
         sl_part.s[i] = -1e10;
   } else {
      sl_part.x    =                 sn / sd   * off_scl[SCLS] + off_scl[OFFS];
      sl_part.y    =                 ln / ld   * off_scl[SCLL] + off_scl[OFFL];
      sl_part.s[2] = (sn*sddx + sd*sndx)/sd/sd * off_scl[SCLS] / off_scl[OFFX];
      sl_part.s[3] = (sn*sddy + sd*sndy)/sd/sd * off_scl[SCLS] / off_scl[OFFY];
      sl_part.s[4] = (sn*sddz + sd*sndz)/sd/sd * off_scl[SCLS] / off_scl[OFFZ];
      sl_part.s[5] = (ln*lddx + ld*lndx)/ld/ld * off_scl[SCLL] / off_scl[OFFX];
      sl_part.s[6] = (ln*lddy + ld*lndy)/ld/ld * off_scl[SCLL] / off_scl[OFFY];
      sl_part.s[7] = (ln*lddz + ld*lndz)/ld/ld * off_scl[SCLL] / off_scl[OFFZ];
   }
}

