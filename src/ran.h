//ran.h

#ifndef __ran__
#define __ran__

#include<cmath>

/*!
  \brief Random number generator

  From Numerical Recipes 3
*/

struct ran {
  unsigned long long int u,v,w;
  ran(unsigned long long int seed=10214L) { set_seed(seed); }

  void set_seed(unsigned long long int seed) {
    v = 4101842887655102017LL;
    w = 1LL;
    u = seed ^ v; int64();
    v = u; int64();
    w = v; int64();
  }

  inline unsigned long long int int64() {
    u = u * 2862933555777941757LL + 7046029254386353087LL;
    v ^= v >> 17; 
    v ^= v << 31; 
    v ^= v >> 8;
    w = 4294957665U*(w & 0xffffffff) + (w >> 32);
    unsigned long long int x = u ^ (u << 21); 
    x ^= x >> 35; 
    x ^= x << 4;
    return (x + v) ^ w;
  }
  inline double doub() { return 5.42101086242752217E-20 * int64(); }
  inline unsigned int int32() { return (unsigned int)int64(); }

  double gauss() {
    double u1,v1,x,y,q;
    do {
      u1 = doub();
      v1 = 1.7156*(doub()-0.5);
      x = u1 - 0.449871;
      y = fabs(v1) + 0.386595;
      q = x*x + y*(0.19600*y - 0.25472*x);
    } while (q > 0.27597 && (q > 0.27846 || v1*v1 > -4.*log(u1)*u1*u1) );
    return v1/u1;
  }
  
};


#endif
