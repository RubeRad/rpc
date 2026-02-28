#ifndef RPC_RPC_H
#define RPC_RPC_H

#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>

using std::string;
using std::istream;
using std::ifstream;
using std::cout;
using std::endl;
using std::vector;

namespace RPC_NS {

typedef int errorType;
#define IS_OK(e)      ((e)>=0)
#define IS_SUCCESS(e) ((e)==0)
#define IS_FAILURE(e) ((e)<0)

enum off_scl_enum {
   OFFX=0,
   OFFY,
   OFFZ,
   OFFS,
   OFFL,
   SCLX,
   SCLY,
   SCLZ,
   SCLS,
   SCLL,
   NUM_OFF_SCL_ENUM
};

// RPC00B coefficient order
enum coeff_enum {
   COEFF_1=0,  // 0 
   COEFF_X,    // 1 
   COEFF_Y,    // 2 
   COEFF_Z,    // 3 
   COEFF_XY,   // 4   RPC00A differences
   COEFF_XZ,   // 5    |
   COEFF_YZ,   // 6    v
   COEFF_XX,   // 7   XYZ  
   COEFF_YY,   // 8   XX
   COEFF_ZZ,   // 9   YY
   COEFF_XYZ,  // 10  ZZ
   COEFF_XXX,  // 11
   COEFF_XYY,  // 12  XXY
   COEFF_XZZ,  // 13  XXZ
   COEFF_XXY,  // 14  XYY
   COEFF_YYY,  // 15
   COEFF_YZZ,  // 16  YYZ
   COEFF_XXZ,  // 17  XZZ
   COEFF_YYZ,  // 18  YZZ
   COEFF_ZZZ,  // 19
   NUM_COEFF_ENUM
};

   
// start with functions using fundamental types that can be turned into openCL kernels
template<typename T, typename GT, typename IT>
void
llh2sl(T* off_scl, // NUM_OFF_SCL_ENUM=10
       T* coeffs,  // NUM_COEFF_ENUM=80
       GT  lon,
       GT  lat,
       GT  hae,
       IT& samp,
       IT& line)
{
   // if GT is double and T is float, this will normalize in double first, and
   // then truncate to float only at the end
   T x( (lon - off_scl[OFFX]) / off_scl[SCLX] );
   T y( (lat - off_scl[OFFY]) / off_scl[SCLY] );
   T z( (hae - off_scl[OFFZ]) / off_scl[SCLZ] );
   T xx=x*x, yy=y*y, zz=z*z, xy=x*y, xz=x*z, yz=y*z;
   
   T sampn=0, sampd=0, linen=0, lined=0;
   // accumulate in reverse order, adding smallest terms first
   for (int i=19; i>=0; --i) {
      T term;
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

   // for well-fitted RPC reasonably in the image domain, the denominators
   // should never cross 0 (or even go negative)
   if (sampd == 0 || lined == 0) {
      samp = -1e10;
      line = -1e10;
   } else {
      // if T is float and IT is double, do all the arithmetic in double
      samp = sampn;
      samp /= sampd;
      samp *= off_scl[SCLS];
      samp += off_scl[OFFS];

      line = linen;
      line /= lined;
      line *= off_scl[SCLL];
      line += off_scl[OFFL];
   }
}


template<typename T>
void
i2g_dlt(T* coeffs,      // NUM_COEFF_ENUM=80
        T s, T l, T z,  // input normalized image coords and normalized Z
        T& x,           // solve for normalized ground X
        T& y)           //                         and Y
{
   // using just the DLT (linear) coefficients plug in Z, solve rest for X,Y
   //      a + bx + cy (+coeff*z)           g + hx + iy (+coeff*z)
   // l = ------------------------     s = ------------------------
   //      d + ex + fy (+coeff*z)           j + kx + my (+coeff*z)
   //
   T a=coeffs[COEFF_1]    + z*coeffs[COEFF_Z],    b=coeffs[COEFF_X],    c=coeffs[COEFF_Y];
   T d=coeffs[COEFF_1+20] + z*coeffs[COEFF_Z+20], e=coeffs[COEFF_X+20], f=coeffs[COEFF_Y+20];
   T g=coeffs[COEFF_1+40] + z*coeffs[COEFF_Z+40], h=coeffs[COEFF_X+40], i=coeffs[COEFF_Y+40];
   T j=coeffs[COEFF_1+60] + z*coeffs[COEFF_Z+60], k=coeffs[COEFF_X+60], m=coeffs[COEFF_Y+60];

   // switch equations to put sample=imgx first; north-up imgs will be like [1,0;0;-1]
   // sj + skx + smy = g + hx + iy ==> (h-sk)x + (i-sm)y = sj-g;
   // ld + lex + lfy = a + bx + cy ==> (b-le)x + (c-lf)y = ld-a;
   T mat00 = h-s*k,   mat01 = i-s*m,   rhs0 = s*j-g;
   T mat10 = b-l*e,   mat11 = c-l*f,   rhs1 = l*d-a;
   T det = mat00*mat11 - mat01*mat10; // north-up will be like -1
   if (det == 0.0) { // Should not happen with well-fitted RPC
      x = -1e-10;
      y = -1e-10;
   } else {
      // 2x2 inverse of matIJ
      T swap = mat00;  // swap the diagonals
      mat00 = mat11;
      mat11 = swap;
      mat10 *= -1;    // negate the off-diagonals
      mat01 *= -1;
      x = (mat00*rhs0 + mat01*rhs1) / det;
      y = (mat01*rhs0 + mat11*rhs1) / det;
#if 0
      // double-check round-trip back
      T testl = (a + b*x + c*y) / (d + e*x + f*y); 
      T tests = (g + h*x + i*y) / (j + k*x + m*y);
      testl -= l;
      tests -= s;
      int stophere=1;
      stophere++;
#endif
   }
}


// hopefully can just get rid of this by using JET/automatic derivative types in
// the above templates
template<typename T>
void
g2ipartials(T* coeffs,     // NUM_COEFF_ENUM=80
            T x, T y, T z, // normalized ground xyz
            T* sl_part)    // *8: x/samp, y/line, then 3+3 partials (normalized)
{
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
       case COEFF_X:   t=x;    tdx=1;    tdy=0;    tdz=0;    break;
       case COEFF_Y:   t=y;    tdx=0;    tdy=1;    tdz=0;    break;
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
      sl_part.x    =                 sn / sd;
      sl_part.y    =                 ln / ld;
      sl_part.s[2] = (sn*sddx + sd*sndx)/sd/sd;
      sl_part.s[3] = (sn*sddy + sd*sndy)/sd/sd;
      sl_part.s[4] = (sn*sddz + sd*sndz)/sd/sd;
      sl_part.s[5] = (ln*lddx + ld*lndx)/ld/ld;
      sl_part.s[6] = (ln*lddy + ld*lndy)/ld/ld;
      sl_part.s[7] = (ln*lddz + ld*lndz)/ld/ld;
   }
}

   

// We could differentiate the type of the normalizers vs the coefficients, but
// practically the normalizers are fairly round values so not much point to
// differentiate. We'll let them be either float or double together
template<typename T>
class RPC {
 public:
    // offsets: lon,lat,hae,samp,line; then scales, same order
   T off_scl[NUM_OFF_SCL_ENUM];
   // same order as WV RPB file, RPC00B TRE
   T coeffs[NUM_COEFF_ENUM*4]; // 20*4=80
   
   RPC() { ; }
  ~RPC() { ; }

   errorType
   init(std::istream& istr)
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
   init(const std::string&  fname)
   {
      //cout << "fname is " << fname << endl;
      ifstream istr(fname.c_str());
      return init(istr);
   }

   template<typename GT, typename IT>
   void
   g2i(GT& gp, IT& ip)
   {
      llh2sl(off_scl, coeffs,
             gp[0], gp[1], gp[2],
             ip[0], ip[1]);
   }
}; // class RPC

}; // namespace RPC
   

// for ground points we will always use double precision
class ground_coord_type {
 public:
   double x,y,z;
   ground_coord_type(double xx=0, double yy=0, double zz=0)
      : x(xx), y(yy), z(zz) { ; }
   double& operator[](int i) {
      switch(i) { case 0: return x;
                  case 1: return y;
                  case 2: return z; }
      return x; // shouldn't ever happen
   }
};

// for image points, maybe double, maybe float.  double is of course sufficient
// precision, but float should be 7 decimal digits precise almost always, which
// is right on the edge of sufficient for satellite imagery, with precision like
// 99999.01 is good, but 100000.1 marginal
template <typename T>
class image_coord_type {
 public:
   T x,y;
   image_coord_type(T xx=0, T yy=0) : x(xx), y(yy) { ; }
   T& operator[](int i) { return (i ? y : x); }
};
typedef image_coord_type<double> imaged_coord_type;
typedef image_coord_type<float>  imagef_coord_type;

#endif
