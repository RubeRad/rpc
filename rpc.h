#ifndef RPC_RPC_H
#define RPC_RPC_H

#include <istream>
#include <string>
#include <vector>

namespace RPC {

typedef int errorType;
#define IS_OK(e)      ((e)>=0)
#define IS_SUCCESS(e) ((e)==0)
#define IS_FAILURE(e) ((e)<0)

typedef double normalizer_type;
typedef double coefficient_type; // TBD sometimes float
typedef double ground_coord_type;
typedef double image_coord_type; // TBD sometimes float

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
          ground_coord_type* llh,  // 3n array of lon,lat,hae (deg/m)
          image_coord_type*  sl);  // preallocated 2n array for samp,line

   // llh is input/output: preallocated 3n, prepopulated with HAEs
   errorType
   slh2ll(size_t  n,
          image_coord_type*  sl,   // 2n array of samp,line
          ground_coord_type* llh); // preallocated 3n array for lon,lat, hae (deg)

}; // class RPC::RPC


template <typename T>
T* first(std::vector<T>& v) { return (&(v[0])); }



void
llh2sl_single(normalizer_type*  off_scl,
              coefficient_type* coeffs,
              ground_coord_type lon,
              ground_coord_type lat,
              ground_coord_type hae,
              image_coord_type& samp,
              image_coord_type& line);
                    

}; // namespace RPC

#endif
