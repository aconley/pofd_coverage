# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(pofd_coverage, 0.3.3, alexander.conley@colorado.edu)
AC_LANG_CPLUSPLUS
AC_CONFIG_SRCDIR([src])
AM_INIT_AUTOMAKE
AC_CONFIG_MACRO_DIR([m4])

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_LIBTOOL
AC_PROG_INSTALL

dnl Initialize CXXFLAGS and CPPFLAGS
dnl CPPFLAGS are preprocessor flags, CXX are compiler flags
dnl Unless --enable-debug is set, then -O3 is used
AC_ARG_WITH(cppflags,
            [  --with-cppflags=VAL     CPPFLAGS value (def: )],
            [CPPFLAGS="$withval"],[CPPFLAGS=""] )
AC_ARG_WITH(cxxflags,[AC_HELP_STRING([--with-cxxflags=VAL],
	[CXXFLAGS value (def: -Wall)])], [ CXXFLAGS=$withval $CXXFLAGS ], 
	[CXXFLAGS="-Wall $CXXFLAGS"])
AC_ARG_ENABLE(debug,
	      AC_HELP_STRING([--enable-debug],
	                     [Add -g to CXXFLAGS]),
              [ CXXFLAGS="-g -O0 $CXXFLAGS"], [ CXXFLAGS="-O3 $CXXFLAGS"])
AC_ARG_ENABLE(profile,
              [  --enable-profile        Add gprof-friendly stuff to CXXFLAGS],
              [ CXXFLAGS="$CXXFLAGS -pg -fprofile-arcs" ],)
AC_ARG_ENABLE(timing,
	      [  --enable-timing         Output timing information],
	      [ AC_DEFINE( TIMING, 1 ) ],)
AC_ARG_ENABLE(openmp,
              [  --enable-openmp         Enable OpenMP],
              [ CXXFLAGS="$CXXFLAGS -fopenmp" ],)
AC_ARG_ENABLE(analyze,
              [  --enable-analyze        Run static checker (if using clang)],
              [ CXXFLAGS="$CXXFLAGS --analyze" ],)

dnl FFTW
AC_LANG_C
AC_ARG_WITH(fftwinc,
	[AC_HELP_STRING([--with-fftwinc=VAL],
	                [VAL=-IVAL, has fftw3 include files ]) ], 
        [fftwinc="$withval/" CPPFLAGS="$CPPFLAGS -I$fftwinc"],
        [fftwinc=] )
AC_ARG_WITH(fftwlib, 
            [  --with-fftwlib=DIR      DIR=-LDIR, has libfftw3 ],
            fftwlib="$withval/",fftwlib= )

dnl FFTW Header
AC_MSG_RESULT( Checking FFTW3 Headers )
AC_CHECK_HEADERS([fftw3.h],,
  [ if test $fftwinc; then
   AC_MSG_ERROR("fftw3 headers not in location specified by --with-fftwinc")
   elif test $fftwlib; then
    AC_CHECK_HEADER($fftwlib/../include/fftw3.h,
     [fftwinc="$fftwlib/../include/" CPPFLAGS="$CPPFLAGS -I$fftwinc"],
      AC_MSG_ERROR("Can't find fftw3 includes in user specified location; use --with-fftwinc"))
    else
      AC_MSG_ERROR("Can't find fftw3 includes; use --with-fftwinc")
   fi ] )

dnl FFTW lib
if test $fftwlib; then
  AC_MSG_RESULT( Checking User specified FFTW3 )
  LDFLAGS="$LDFLAGS -L$fftwlib"
  LIBS="$LIBS -lfftw3"
  AC_TRY_LINK( [#include<fftw3.h>],
  	       [fftw_plan plan;],
	       [echo "Found libFFTW3"],
	       [AC_MSG_ERROR("FFTW3 library not found in user specified location")])
dnl  AC_CHECK_LIB( fftw3, main, ,
dnl      AC_MSG_ERROR("libfftw3 not in location specified by --with-fftwlib: $fftwlib"),
dnl      [ -lfftw3 -L$fftwlib ] )
else
  AC_MSG_RESULT( Checking FFTW3 )
  AC_CHECK_LIB( fftw3, main, ,
      AC_MSG_ERROR("libfftw3 not found; use --with-fftwlib"), [-lm -lfftw3] )
fi

dnl GSL
AC_LANG_C
AC_ARG_WITH(gsllib, 
            [  --with-gsllib=DIR       DIR=-LDIR, has libgsl and libgslcblas],
            gsllib="$withval/",gsllib= )
if test $gsllib; then
  AC_MSG_RESULT( Checking User specified GSL )
  LDFLAGS="$LDFLAGS -L$gsllib"
  AC_CHECK_LIB( gsl, main, ,
      AC_MSG_ERROR("libgsl not in location specified by --with-gsllib: $gsllib"),
      [ -lm -lgslcblas -L$gsllib ] )
  AC_CHECK_LIB( gslcblas, main, ,
      AC_MSG_ERROR("libgslcblas not in location specified by --with-gsllib: $gsllib"),
      [ -lm -L$gsllib ] )
else
  AC_MSG_RESULT( Checking GSL )
  AC_CHECK_LIB( gsl, main, ,
      AC_MSG_ERROR("libgsl not found; use --with-gsllib"), [-lm -lgslcblas] )
  AC_CHECK_LIB( gslcblas, main, ,
      AC_MSG_ERROR("libgslcblas not found; use --with-gsllib"), [-lm] )
fi

AC_ARG_WITH(gslinc,
        [  --with-gslinc=VAL       VAL=-IVAL, has gsl include files ],
	[gslinc="$withval/" CPPFLAGS="$CPPFLAGS -I$gslinc"],
        [gslinc=] )
AC_MSG_RESULT( Checking GSL Headers )
AC_CHECK_HEADERS([gsl/gsl_errno.h gsl/gsl_math.h gsl/gsl_min.h ],,
  [ if test $gslinc; then
   AC_MSG_ERROR("gsl headers not in location specified by --with-gslinc")
   elif test $gsllib; then
    AC_CHECK_HEADER($gsllib/../include/gsl/gsl_math.h,
     [gslinc="$gsllib/../include/" CPPFLAGS="$CPPFLAGS -I$gslinc"],
      AC_MSG_ERROR("Can't find gsl includes in user specified location; use --with-gslinc"))
    else
      AC_MSG_ERROR("Can't find gsl includes; use --with-gslinc")
   fi ] )


dnl Handle cfitsio library stuff
AC_LANG_C
dnl We'll need to possible handle a user specified CFITSIO location
dnl We can't use --with-cfitsio because autoconf gets confused
dnl with the library name for some strange reason
AC_ARG_WITH( cfitsiolib, 
        [  --with-cfitsiolib=DIR   DIR=-LDIR, has libcfitsio.a],
        cfitsiolib="$withval/" , cfitsiolib= )
if test $cfitsiolib; then
   LDFLAGS="$LDFLAGS -L$cfitsiolib"

fi
AC_CHECK_LIB(cfitsio,main,,
         [AC_MSG_ERROR("Can't find libcfitsio; use --with-cfitsiolib")], 
	    [-lm "-L$cfitsiolib"])
dnl And cfitsio headers
AC_ARG_WITH(cfitsioinc,
     [  --with-cfitsioinc=VAL   VAL=-IVAL, has cfitsio include files ],
	[cfitsioinc="$withval/" CPPFLAGS="$CPPFLAGS -I$cfitsioinc"],
     [cfitsioinc=] )
AC_CHECK_HEADER(fitsio.h,,
  [ if test $cfitsioinc; then
   AC_MSG_ERROR("cfitsio headers not in location specified by --with-cfitsioinc")
   elif test $cfitsiolib; then
    AC_CHECK_HEADER($cfitsiolib/../include/fitsio.h,
     [cfitsioinc="$cfitsiolib/../include/" CPPFLAGS="$CPPFLAGS -I$cfitsioinc"],
      AC_MSG_ERROR("Can't find cfitsio includes; use --with-cfitsioinc"))
    else
      AC_MSG_ERROR("Can't find cfitsio includes; use --with-cfitsioinc")
   fi ] )

dnl Tests
AC_LANG_CPLUSPLUS
AM_CONDITIONAL(DOTEST, false)
AC_ARG_ENABLE(test,
        [  --enable-test=DIR       Enable testing],
        usetest=$enableval, usetest="")
if test $usetest; then
   AC_MSG_RESULT(Enabling testing, unit and other)
   AM_CONDITIONAL(DOTEST, true)

   dnl Find gtest library
   AC_ARG_WITH(gtestlib,
               [  --with-gtestlib=DIR       DIR=-LDIR, has libgtest],
               gtestlib="$withval/",gtestlib= )
   if test $gtestlib; then
      AC_MSG_RESULT( Checking User Specified gtest )
      LDFLAGS="$LDFLAGS -L$gtestlib"
      AC_CHECK_LIB( gtest_main, testing::InitGoogleTest, ,
                    AC_MSG_ERROR("libgtest not in location specified by --with-gtestlib: $gtestlib"),
                    [ -lgtest -lgtest_main ] )
   else
      AC_MSG_RESULT( Checking gtest_main )
      AC_CHECK_LIB( gtest_main, main, ,
      AC_MSG_ERROR("libgtest_main not found -- use --with-gtestlib"),
                   [ -lgtest -lgtest_main ] )
      AC_MSG_RESULT( Checking gtest )
      AC_CHECK_LIB( gtest, main, ,
      AC_MSG_ERROR("libgtest not found -- use --with-gtestlib"),
                   [ -lgtest -lgtest_main ] )
   fi

   dnl Find gtest includes
   AC_ARG_WITH(gtestinc,
        [  --with-gtestinc=VAL      VAL=-IVAL, has gtest include files ],
        [gtestinc="$withval/" CPPFLAGS="$CPPFLAGS -I$gtestinc"],
        [gtestinc=] )
   AC_MSG_RESULT( Checking gtest Headers )
   AC_CHECK_HEADERS([gtest/gtest.h],,
     [ if test $gtestinc; then
         AC_MSG_ERROR("gtest headers not in location specified by --with-gtestinc")
       elif test $gtestlib; then
         AC_CHECK_HEADER($gtestlib/../include/test.h,
            [gtestinc="$gtestlib/../include/" CPPFLAGS="$CPPFLAGS -I$gtestinc"],
             AC_MSG_ERROR("Can't find gtest includes in user specified location; use --with-gtestinc"))
       else
          AC_MSG_ERROR("Can't find gtest includes; use --with-gtestinc")
       fi ] )
fi

dnl HDF5
AC_LANG_C
AC_ARG_WITH(hdfinc,
        [AC_HELP_STRING([--with-hdfinc=VAL],
                        [VAL=-IVAL, has HDF5 include files ]) ], 
        [hdfinc="$withval/" CPPFLAGS="$CPPFLAGS -I$hdfinc"],
        [hdfinc=] )
AC_ARG_WITH(hdflib,
	[AC_HELP_STRING([--with-hdflib=DIR],
			[DIR=-LDIR, has libhdf5])],
	hdflib="$withval/", hdflib= )

dnl HDF5 Header
AC_MSG_RESULT( Checking HDF5 Headers )
AC_CHECK_HEADERS([hdf5.h],,
  [ if test $hdfinc; then
   AC_MSG_ERROR("hdf5 headers not in location specified by --with-hdfinc")
   elif test $hdflib; then
    AC_CHECK_HEADER($hdflib/../include/hdf5.h,
     [hdfinc="$hdflib/../include/" CPPFLAGS="$CPPFLAGS -I$hdfinc"],
      AC_MSG_ERROR("Can't find HDF5 includes in user specified location; use --with-hdfinc"))
    else
      AC_MSG_ERROR("Can't find HDF5 includes; use --with-hdf5inc")
   fi ] )
dnl HDF5 lib
if test $hdflib; then
  AC_MSG_RESULT( Checking User specified HDF5 )
  LDFLAGS="$LDFLAGS -L$hdflib"
  LIBS="$LIBS -lhdf5"
  AC_CHECK_LIB( hdf5, main, ,
      AC_MSG_ERROR("libhdf5 not in location specified by --with-hdflib: $hdflib"),
      [ -lhdf5 -L$hdflib ] )

else
  AC_MSG_RESULT( Checking HDF5 )
  AC_CHECK_LIB( hdf5, main, ,
      AC_MSG_ERROR("libhdf5 not found; use --with-hdflib"), )
fi

# Checks for header files.
AC_HEADER_STDC

AC_CONFIG_FILES([Makefile src/Makefile test/Makefile])
AC_OUTPUT
