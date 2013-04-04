#include<iostream>
#include<cmath>
#include<limits>
#include<cstdlib>
#include<algorithm>
#include<sstream>

#include "../include/utility.h"
#include "../include/pofdExcept.h"

/*! 
  Breaks an input string up into a vector of string, which correspond
  to the input string split on spaces.  Cheerfully stolen from Rob
  Knop.
*/
void utility::stringwords(const std::string &ins,
			  std::vector<std::string> &words) {
  std::string s,tmp;
  unsigned int i,p;
  int first,last;

  s = ins;

  // Trim spaces from beginning and end

  first=s.find_first_not_of(" ");
  if (first==-1) {
    s="";
  } else {
    last=s.find_last_not_of(" ");
    s=s.substr(first,last-first+1);
  }

  words.clear();

  p=s.find_first_not_of(" \t\r\n");
  if (p>=s.size()) return;
  while ((i=s.find_first_of(" \t\r\n",p))) {
    tmp=s.substr(p,i-p);
    words.push_back(tmp);
    p=s.find_first_not_of(" \t\r\n",i);
    if (p>=s.size()) return;
  }
  tmp=s.substr(p);
  words.push_back(tmp);
}

/*!
  Finds the position in a sorted array of the last element
  less than a specified value. ndata is returned if a problem is
  encountered.
 */
unsigned int utility::binary_search_lt(double value, double* data, 
				       unsigned int ndata) {
  unsigned int l, u; //Lower and upper bounds of current box
  unsigned int m; //Current element

  const unsigned int maxrep = 100; //Assume we will never have 2^100 elements

  //Quick exit checks
  if (ndata <= 1) return 0;
  if (value > data[ndata-1]) return ndata-1;
  if (value <= data[0]) return ndata;

  l = 0; u = ndata-1;
  if (data[l] > data[u]) return ndata;
  
  //Main loop
  unsigned int i;
  for (i = 0; i < maxrep; ++i) {
    m = (l+u)/2;
    if (l == m) break; //l is what we want.
    if (data[m] < value) l=m; else u=m;
  }
  if (i == maxrep) return ndata;
  return l;
}

/*!
  Finds the position in a sorted array of the last element
  less than a specified value. ndata is returned if a problem is
  encountered.
 */
unsigned int utility::binary_search_lt(double value, 
				       const std::vector<double>& data) {
  unsigned int l, u; //Lower and upper bounds of current box
  unsigned int m; //Current element

  const unsigned int maxrep = 100; //Assume we will never have 2^100 elements

  //Quick exit checks
  unsigned int ndata = data.size();
  if (ndata <= 1) return 0;
  if (value > data[ndata-1]) return ndata-1;
  if (value <= data[0]) return ndata;

  l = 0; u = ndata-1;
  if (data[l] > data[u]) return ndata;
  
  //Main loop
  unsigned int i;
  for (i = 0; i < maxrep; ++i) {
    m = (l+u)/2;
    if (l == m) break; //l is what we want.
    if (data[m] < value) l=m; else u=m;
  }
  if (i == maxrep) return ndata;
  return l;
}

unsigned int utility::binary_search_lte(double value, double* data,
					unsigned int ndata) {
  unsigned int idx = binary_search_gt(value,data,ndata);
  if (idx == 0) return ndata; else return idx-1;
}

/*!
  Finds the position in a sorted array of the first element
  greater than a specified value. ndata is returned if a problem is
  encountered.
 */
unsigned int utility::binary_search_gt(double value, double* data, 
				       unsigned int ndata) {
  unsigned int l, u; //Lower and upper bounds of current box
  unsigned int m; //Current element

  const unsigned int maxrep = 100; //Assume we will never have 2^100 elements

  //Quick exit checks
  if (ndata <= 0) return ndata;
  if (ndata == 1) {
    if (data[0] > value) return 0; else return ndata;
  }
  if (value >= data[ndata-1]) return ndata;
  if (value < data[0]) return 0;

  l = 0; u = ndata-1;
  if (data[l] > data[u]) return ndata;
  
  //Main loop
  unsigned int i;
  for (i = 0; i < maxrep; ++i) {
    m = (l+u)/2;
    if (l == m) break; //u is what we want.
    if (data[m] > value) u=m; else l=m;
  }
  if (i == maxrep) return ndata;
  return u;
}

/*!
  Finds the position in a sorted array of the first element
  greater than a specified value. ndata is returned if a problem is
  encountered.
*/
  
unsigned int utility::binary_search_gt(double value,std::vector<double> data){
  unsigned int l, u; //Lower and upper bounds of current box
  unsigned int m; //Current element
  unsigned int ndata, i; //Length of data array
  const unsigned int maxrep = 100; //Assume we will never have 2^100 elements
  ndata = data.size();
  if (ndata == 1) {
    if (data[0] > value) return 0; else return ndata;
  }
  if (value >= data[ndata-1]) return ndata;
  if (value < data[0]) return 0;
  l = 0; u = ndata-1;
  if (data[l] > data[u]) return ndata;
  for (i = 0; i < maxrep; ++i) {
    m = (l+u)/2;
    if (l == m) break; //u is what we want.
    if (data[m] > value) u=m; else l=m;
  }
  if (i == maxrep) return ndata;
  return u;
}

/*!
  Finds the position in a reverse sorted array of the last element
  greater than a specified value. ndata is returned if a problem is
  encountered.
 */
int utility::binary_search_rev(double value, double* data, 
			       unsigned int ndata) {
  unsigned int l, u; //Lower and upper bounds of current box
  unsigned int m; //Current element

  const unsigned int maxrep = 100; //Assume we will never have 2^100 elements

  //Quick exit checks
  if (ndata <= 1) return ndata;
  if (value < data[ndata-1]) return ndata-1;
  if (value < data[0]) return ndata;

  l = 0; u = ndata-1;
  if (data[l] < data[u]) return -1;
  
  //Main loop
  unsigned int i;
  for (i = 0; i < maxrep; ++i) {
    m = (l+u)/2;
    if (l == m) break; //l is what we want.
    if (data[m] > value) l=m; else u=m;
  }
  if (i == maxrep) return ndata;
  return l;
}

/*!
  Not for the faint of heart, and assumes unsigned int is
  32 bits.

  From http://graphics.stanford.edu/~seander/bithacks.html#IntegerLog
 */
unsigned int utility::log2( unsigned int val ) {
  if (val == 1) {
    return 0;
  } else if ( val & (val-1) ) {
    //Not a power of 2
    const unsigned int b[] = {0x2, 0xC, 0xF0, 0xFF00, 0xFFFF0000};
    const unsigned int S[] = {1, 2, 4, 8, 16};

    register unsigned int r = 0;
    for (int i = 4; i >= 0; i--) {
      if (val & b[i]) {
	val >>= S[i];
	r |= S[i];
      } 
    }
    return r;
  } else {
    //A power of 2
    const unsigned int b[] = {0xAAAAAAAA, 0xCCCCCCCC, 0xF0F0F0F0, 
			      0xFF00FF00, 0xFFFF0000};
    unsigned int r = (val & b[0]) != 0;
    for (int i = 4; i > 0; i--)
      r |= ((val & b[i]) != 0) << i;
    return r;
  }
}

double utility::logfactorial(double value) {
  //Lanzcos (1964), SIAM Journal on Numerical Analysis, ser B, vol 1., p 86
  if (value < 0) return std::numeric_limits<double>::quiet_NaN();
  if (value == 0 || value == 1) return 0.0;
  const unsigned int nterms = 14;
  double vp1, temp, y, series;
  const double coeffs[nterms] = {57.1562356658629235,
				 -59.5979603554754915,
				 14.1360979747417471,
				 -0.491913816097620199,
				 0.339946499848118887e-4,
				 0.465236289270485756e-4,
				 -0.983744753048795646e-4,
				 0.158088703224912494e-3,
				 -0.210264441724104883e-3,
				 0.217439618115212643e-3,
				 -0.164318106536763890e-3,
				 0.84412239838527433e-4,
				 -0.261908384015814087e-4,
				 0.368991826595316234e-5};
  vp1 = value+1.0;
  temp = vp1 + 5.2421875;
  temp = (vp1+0.5)*log(temp)-temp;
  series = 0.999999999999997092;
  y = vp1;
  for (unsigned int i = 0; i < nterms; ++i) series += coeffs[i]/++y;
  return temp + log(2.5066282746310005*series/vp1);
}
