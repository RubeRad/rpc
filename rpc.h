#ifndef RPC_RPC_H
#define RPC_RPC_H

#include <istream>
#include <string>
#include <vector>

#include <CL/cl.h>

namespace RPC {

typedef int errorType;
#define IS_OK(e)      ((e)>=0)
#define IS_SUCCESS(e) ((e)==0)
#define IS_FAILURE(e) ((e)<0)

typedef double  normalizer_type;
typedef double  coefficient_type; // TBD sometimes float
typedef cl_double3 ground_coord_type;
typedef cl_double2 image_coord_type; // TBD sometimes float


enum scl_off_enum {
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
};
enum coeff_enum {
   COEFF_1=0,  // 0
   COEFF_X,    // 1
   COEFF_Y,    // 2
   COEFF_Z,    // 3
   COEFF_XY,   // 4
   COEFF_XZ,   // 5
   COEFF_YZ,   // 6
   COEFF_XX,   // 7
   COEFF_YY,   // 8
   COEFF_ZZ,   // 9
   COEFF_XYZ,  // 10
   COEFF_XXX,  // 11
   COEFF_XYY,  // 12
   COEFF_XZZ,  // 13
   COEFF_XXY,  // 14
   COEFF_YYY,  // 15
   COEFF_YZZ,  // 16
   COEFF_XXZ,  // 17
   COEFF_YYZ,  // 18
   COEFF_ZZZ,  // 19
};
   

class RPC {
 public:
   normalizer_type  off_scl[10]; // offsets: lon,lat,hae,samp,line
                                 // then scales, same order
   coefficient_type coeffs[80];  // same order as WV RPB file
   
   RPC() { ; }
  ~RPC() { ; }

   errorType init(      std::istream& istr);
   errorType init(const std::string&  fname);

   errorType
   llh2sl(size_t  n,
          ground_coord_type* llh,  // n-array of lon,lat,hae (deg/m) triplets
          image_coord_type*  sl)   // preallocated n-array for samp,line duplets
   const;

   // llh is input/output: preallocated, prepopulated with HAEs
   errorType
   slh2ll(size_t  n,
          image_coord_type*  sl,   // n-array of samp,line duplets
          ground_coord_type* llh)  // preallocated n-array for lon,lat,hae (deg)
   const;

}; // class RPC::RPC


template <typename T>
T* first(std::vector<T>& v) { return (&(v[0])); }



void
llh2sl_single(const normalizer_type*  off_scl,
              const coefficient_type* coeffs,
              double  lon,
              double  lat,
              double  hae,
              double& samp,
              double& line);

void
g2ipartials(normalizer_type* off_scl,
            coefficient_type* coeffs,
            const cl_double3& llh,
            cl_double8& sl_part);

void
i2g_dlt(const coefficient_type* coeffs,
        double nix, double niy, // normalized image coords
        double ngz,             // normalized ground Z
        double& ngx,  // solve for normalized ground X
        double& ngy); //                         and Y
                    

}; // namespace RPC

#endif
