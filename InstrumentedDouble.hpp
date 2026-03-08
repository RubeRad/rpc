#ifndef INSTRUMENTED_DOUBLE_HPP
#define INSTRUMENTED_DOUBLE_HPP

class InstrumentedDouble {
 public:
   double x;

   static long PLUS_COUNT; // and minus
   static long MULT_COUNT; // and divide
   static void reset() { PLUS_COUNT = MULT_COUNT = 0; }

   InstrumentedDouble(double y=0.0) : x(y)   { ; }
   bool operator==(const InstrumentedDouble& y) { return (x == y.x); }
   bool operator< (const InstrumentedDouble& y) { return (x < y.x); }

   InstrumentedDouble   operator+ (const InstrumentedDouble& y) const { PLUS_COUNT++; return InstrumentedDouble(x+y.x); }
   InstrumentedDouble   operator- (const InstrumentedDouble& y) const { PLUS_COUNT++; return InstrumentedDouble(x-y.x); }
   InstrumentedDouble   operator* (const InstrumentedDouble& y) const { MULT_COUNT++; return InstrumentedDouble(x*y.x); }
   InstrumentedDouble   operator/ (const InstrumentedDouble& y) const { MULT_COUNT++; return InstrumentedDouble(x/y.x); }
   InstrumentedDouble&  operator+=(const InstrumentedDouble& y)       { PLUS_COUNT++; x+=y.x; return *this; }
   InstrumentedDouble&  operator-=(const InstrumentedDouble& y)       { PLUS_COUNT++; x-=y.x; return *this; }
   InstrumentedDouble&  operator*=(const InstrumentedDouble& y)       { MULT_COUNT++; x*=y.x; return *this; }
   InstrumentedDouble&  operator/=(const InstrumentedDouble& y)       { MULT_COUNT++; x/=y.x; return *this; }
};

// ISO C++ forbids in-class initialization of non-const static member
long InstrumentedDouble::PLUS_COUNT=0;
long InstrumentedDouble::MULT_COUNT=0;

#endif
