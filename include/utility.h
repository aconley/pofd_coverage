//utility.h

#ifndef __utility__
#define __utility__

#include<vector>
#include<string>

/*!
  \brief Utility functions
*/
namespace utility {

  /* \brief Break input string up into words */
  void stringwords(const std::string &ins,
		   std::vector<std::string> &words);

  /* \brief Binary search for last value less than specified amount */
  unsigned int binary_search_lt(double value, double* data, unsigned int ndata);

  /* \brief Binary search for last value less than specified amount */
  unsigned int binary_search_lt(double value, const std::vector<double>&);

  /* \brief Binary search for first value greater than specified amount */
  unsigned int binary_search_gt(double value, double* data, unsigned int ndata);

  /*! \brief Binary search for first value greater than specified amount */
  unsigned int binary_search_gt(double value,std::vector<double> data);

  /* \brief Binary search for last value less than or equal to 
     specified amount */
  unsigned int binary_search_lte(double value, double* data, 
				 unsigned int ndata);

  /* \brief Binary search for last value greater than specified amount in reverse array*/
  int binary_search_rev(double value, double* data, 
			unsigned int ndata);

  /* \brief Binary logarithm (\f$log_{2}\f$) of 32 bit integer */
  unsigned int log2( unsigned int );

  /* \brief Log of factorial of argument */
  double logfactorial( double );

  std::vector<double> translateParams(const std::vector<double>& knots, 
				      const std::vector<double>& knotvals);

  /* \brief Abstract base interpolation/extrapolation class */
  class base_interp {
  protected :
    double *c;
    double *d;
    unsigned int M;
  public :
    base_interp();
    base_interp(unsigned int );
    ~base_interp();
    unsigned int getM() const {return M;} //!< Get's poly order+1
    double getC(unsigned int i) const { return c[i]; }
    double getD(unsigned int i) const { return d[i]; }
    void resize( unsigned int );
    virtual double interpol(unsigned int, double, unsigned int,
			    double*, double*, double&, bool&) const = 0;
  };

  /* \brief Polynomial interpolation/extrapolation */
  class poly_interp : public base_interp {
  public :
    poly_interp();
    poly_interp( unsigned int );
    ~poly_interp();
    double interpol(unsigned int, double, unsigned int,
		    double*, double*,double&, bool&) const;
  };

  /* \brief Rational interpolation/extrapolation */
  class rat_interp : public base_interp {
  public :
    rat_interp();
    rat_interp( unsigned int );
    ~rat_interp();
    double interpol(unsigned int, double, unsigned int,
		    double*, double*,double&, bool&) const;
  };

}

#endif
