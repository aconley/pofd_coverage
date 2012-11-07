#include<paramSet.h>
#include<pofdExcept.h>

paramSet::paramSet(unsigned int NPARAMS, unsigned int NNOISEPARAMS) {
  nparams = nnoiseparams = 0;
  paramvals = noiseParams = NULL;
  setNParams(NPARAMS);
  setNNoiseParams(NNOISEPARAMS);
}

paramSet::paramSet(const std::vector<double>& vec, 
		   const std::vector<double>& nvec) {
  nparams = nnoiseparams = 0;
  paramvals = noiseParams = NULL;
  setNParams(vec.size());
  setNNoiseParams(nvec.size());
  if (nparams > 0) for (unsigned int i = 0; i < nparams; ++i)
		     paramvals[i]=vec[i];
  if (nnoiseparams > 0) for (unsigned int i = 0; i < nnoiseparams; ++i)
			  noiseParams[i] = nvec[i];
}

paramSet::~paramSet() {
  if (paramvals != NULL) delete[] paramvals;
  if (noiseParams != NULL) delete[] noiseParams;
}

void paramSet::setNParams(unsigned int npar) {
  if (npar == nparams) return;
  if (paramvals != NULL) delete[] paramvals;
  if (npar > 0) paramvals = new double[npar]; else paramvals=NULL;
  nparams=npar;
}

void paramSet::setNNoiseParams(unsigned int npar) {
  if (npar == nnoiseparams) return;
  if (noiseParams != NULL) delete[] noiseParams;
  if (npar > 0) noiseParams = new double[npar]; else noiseParams=NULL;
  nnoiseparams=npar;
}

paramSet& paramSet::operator=(const paramSet& other) {
  if (this == &other) return *this;
  setNParams(other.nparams);
  for (unsigned int i = 0; i < other.nparams; ++i)
    paramvals[i] = other.paramvals[i];
  setNNoiseParams(other.nnoiseparams);
  for (unsigned int i = 0; i < other.nnoiseparams; ++i)
    noiseParams[i] = other.noiseParams[i];
  return *this;
}

/*!
  \param[in] vec Input parameter vector
  Doesn't allow for resizing, doesn't change noise values
 */
void paramSet::setParamValues(const std::vector<double>& vec) {
  if (vec.size() != nparams)
    throw pofdExcept("paramSet","setParamValues",
		     "Input vector wrong length",1);
  for (unsigned int i = 0; i < nparams; ++i)
    paramvals[i]=vec[i];
}


/*!
  \param[in] vec Input parameter vector
  \param[in] nvec Input noise parameter vector
  Doesn't allow for resizing
 */
void paramSet::setParamValues(const std::vector<double>& vec, 
			      const std::vector<double>& nvec) {
  if (vec.size() != nparams)
    throw pofdExcept("paramSet","setParamValues",
		     "Input vector wrong length",1);
  if (nvec.size() != nnoiseparams)
    throw pofdExcept("paramSet","setParamValues",
		     "Input noise vector wrong length",2);
  for (unsigned int i = 0; i < nparams; ++i)
    paramvals[i]=vec[i];
  for (unsigned int i = 0; i < nparams; ++i)
    noiseParams[i]=nvec[i];
}

/*!
  Doesn't allow resizing
*/
void paramSet::readFromStream(std::istream& is) {
  for (unsigned int i = 0; i < nparams; ++i)
    is >> paramvals[i];
  for (unsigned int i = 0; i < nnoiseparams; ++i)
    is >> noiseParams[i];
}

bool paramSet::writeToStream( std::ostream& os ) const {
  for (unsigned int i = 0; i < nparams; ++i)
    os << "   " << paramvals[i];
  for (unsigned int i = 0; i < nnoiseparams; ++i)
    os << "   " << noiseParams[i];

  return true;
}

std::istream& operator>>(std::istream& is, paramSet& p) {
  p.readFromStream(is);
  return is;
}

std::ostream& operator<<(std::ostream& os, const paramSet& p) {
  p.writeToStream(os);
  return os;
}

