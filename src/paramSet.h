//paramset.h

#ifndef __paramset__
#define __paramset__

#include<vector>
#include<istream>
#include<ostream>

/*!
  \brief Class for holding parameter values of the model, and
  noise description parameters
*/
class paramSet {
 private:
  unsigned int nparams;
  double* paramvals;
  unsigned int nnoiseparams;
  double* noiseParams;
 public:
  paramSet(unsigned int NPARAMS=0, unsigned int NNOISEPARAMS=1);
  paramSet(const std::vector<double>&, const std::vector<double>&);
  ~paramSet();

  void setNParams(unsigned int);
  void setNNoiseParams(unsigned int);

  unsigned int getNParams() const { return nparams; }
  unsigned int getNNoiseParams() const { return nnoiseparams; }
  const double& operator[](unsigned int i) const { return paramvals[i]; }
  double& operator[](unsigned int i) { return paramvals[i]; }

  void setParamValues(const std::vector<double>&, const std::vector<double>&);
  void setParamValues(const std::vector<double>&);
  void setParamValue(unsigned int i, double val) { paramvals[i]=val; }
  void setNoiseParamValue(unsigned int i, double val) { noiseParams[i]=val;}
  paramSet& operator=(const paramSet&);

  double getNoiseParam(unsigned int idx) const { return noiseParams[idx]; }
  void setNoiseParam(unsigned int idx, double val) { noiseParams[idx] = val; }

  //Input
  void readFromStream(std::istream& is);

  //Output
  bool writeToStream(std::ostream& os) const;
  
};

std::istream& operator>>(std::istream& is, paramSet& p);
std::ostream& operator<<(std::ostream& os, const paramSet& p);

#endif
