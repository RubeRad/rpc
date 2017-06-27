#ifndef RPC_RPC_H
#define RPC_RPC_H

#include <istream>
#include <string>

namespace RPC {

typedef int errorType;
#define IS_OK(e)      ((e)>=0)
#define IS_SUCCESS(e) ((e)==0)
#define IS_FAILURE(e) ((e)<0)

typedef double point2[2]; // ALWAYS X,Y   (samp,line)
typedef double point3[3]; // ALWAYS X,Y,Z (lon, lat, hae)
typedef double normalizer_type;
typedef double coefficient_type; // TBD sometimes float

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
          point3* llh,  // lon,lat,hae (deg/m)
          point2* sl);  // samp,line

   // HAE prepopulated in input/output param llh!
   errorType
   slh2ll(size_t  n,
          point2* sl,   // samp,line
          point3* llh); // lon,lat, hae (deg)

}; // class RPC::RPC

}; // namespace RPC

#endif
