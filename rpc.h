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
          image_coord_type*  sl);  // preallocated n-array for samp,line duplets

   // llh is input/output: preallocated, prepopulated with HAEs
   errorType
   slh2ll(size_t  n,
          image_coord_type*  sl,   // n-array of samp,line duplets
          ground_coord_type* llh); // preallocated n-array for lon,lat,hae (deg)

}; // class RPC::RPC


template <typename T>
T* first(std::vector<T>& v) { return (&(v[0])); }



void
llh2sl_single(normalizer_type*  off_scl,
              coefficient_type* coeffs,
              double  lon,
              double  lat,
              double  hae,
              double& samp,
              double& line);
                    

}; // namespace RPC

#endif
