#include<iostream>
#include<sstream>

#include<pofdExcept.h>

//Cheerfully stolen from deepexcpt, written by Rob Knop

/*!
  Does the work of the constructors, except for setting the flags
  \param[in] inclass Class generating exception
  \param[in] inmethod Method generating exception
  \param[in] inerrstr Error string 
  \param[in] inerr Error number
*/
void pofdExcept::init(const std::string& inclass,
		      const std::string& inmethod,
		      const std::string& inerrstr,int inerr) {
  errclass=inclass;
  errmethod=inmethod;
  errstr=inerrstr;
  errnum=inerr;
}


pofdExcept::pofdExcept() {
  init("","","",0);
  classset = methodset = strset = errset = false;
}

/*!
  Most basic error, specifying only the error message
 */
pofdExcept::pofdExcept(const std::string errstr) {
  init("","",errstr,0);
  strset = true;
}

/*!
  Error with error string and error number
*/
pofdExcept::pofdExcept(const std::string errstr,int err) {
  init("","",errstr,err);
  strset = errset = true;
}

/*!
  Error with error string, class and method generating exception
*/
pofdExcept::pofdExcept(const std::string errclass,
				 const std::string errmethod,
				 const std::string errstr) {
  init(errclass,errmethod,errstr,0);
  classset = methodset = strset = true;
}

/*!
  Full error specification: error std::string, number, class, and method.
 */
pofdExcept::pofdExcept(const std::string errclass,
				 const std::string errmethod,
				 const std::string errstr,int err) {
  init(errclass,errmethod,errstr,err);
  classset = methodset = strset = errset = true;
}

std::string pofdExcept::what() const {
  std::stringstream str;
  bool first; //Last element won't have a linebreak at the end
  first = true;
  if (classset) {
    str << "Error Class/Namespace: " << errclass;
    first = false;
  }
  if (methodset) {
    if (!first) str << std::endl;
    str << "Method: " << errmethod;
    first = false;
  }
  if (strset) {
    if (!first) str << std::endl;
    str << "Error Message: " << errstr;
    first = false;
  }
  if (errset) {
    if (!first) str << std::endl;
    str << "Error Code: " << errnum;
  }
  return str.str();
}

/*
  Provides output capabilities, and is smart enough not to output
  information that hasn't been set.
 */
std::ostream& operator<<(std::ostream& os, const pofdExcept& err) {
  os << err.what();
  return os;
}
