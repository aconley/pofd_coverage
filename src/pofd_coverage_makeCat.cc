#include<iostream>
#include<string>
#include<iomanip>
#include<vector>
#include<fstream>
#include<cstdlib>
#include<ctime>
#include<utility>

#include<getopt.h>

#include<fitsio.h>

#include "../include/numberCounts.h"
#include "../include/numberCountsDouble.h"
#include "../include/ran.h"
#include "../include/pofdExcept.h"
#include "../include/global_settings.h"

#pragma GCC diagnostic ignored "-Wwrite-strings"

//Set up global option index that can be used for both single and double case
static struct option long_options[] = {
  {"double",no_argument,0,'d'},
  {"help", no_argument, 0, 'h'},
  {"seed", required_argument, 0,'S'},
  {"verbose",no_argument,0,'v'},
  {"version",no_argument,0,'V'},
  {0,0,0,0}
};
char optstring[] = "dhS:vV";

int makeCatSingle(int argc, char **argv) {

  unsigned int n0;
  std::string modelfile, outputfile; 
  unsigned long long int seed;
  bool verbose, have_user_seed;

  //Defaults
  verbose             = false;
  have_user_seed      = false;
  seed                = 0;


  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'S' :
      have_user_seed = true;
      seed = static_cast<unsigned long long int>( atoi(optarg) );
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc-2 ) {
    std::cout << "Some required arguments missing" << std::endl;
    std::cout << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  modelfile  = std::string(argv[optind]);
  n0         = atoi(argv[optind + 1]);
  outputfile = std::string(argv[optind + 2]);

  if (n0 == 0) {
    std::cout << "Invalid (non-positive) n0: " << n0 << std::endl;
    return 1;
  }

  try {

    //Set up random number generator
    ran rangen;
    if (!have_user_seed) {
      seed = static_cast<unsigned long long int>(time(NULL));
      seed += static_cast<unsigned long long int>(clock());
    }
    rangen.set_seed(seed);

    //Set up the model
    numberCounts model(modelfile);
    if (!model.isValid())
      throw pofdExcept("pofd_coverage_makeCat", "makeCatSingle",
		       "Trying to realize model with invalid parameters", 1);

    //Get the fluxes
    float *flux;
    flux = new float[n0];
    for (unsigned int i = 0; i < n0; ++i)
      flux[i] = model.genSource(rangen.doub());

    if (verbose) std::cout << "Writing catalog to " << outputfile 
			   << std::endl;

    //Do the write
    int status = 0;
    fitsfile *fp;
    fits_create_file(&fp, outputfile.c_str(), &status);
    if (status) {
      fits_report_error(stderr,status);
      return status;
    }

    int tfields = 1;
    long nrows = static_cast<long>(n0);
    char extname[] = "CATALOG";
    //These make g++ complain, but there isn't much I can do about it
    // because I can't change the call signature of fits_create_tbl
    // to make them const*
    char *ttype[] = {"FluxDensity"};
    char *tform[] = {"1E"};
    char *tunit[] = {"Jy"};
    fits_create_tbl(fp, BINARY_TBL, nrows, tfields, ttype, tform,
		    tunit, extname, &status);
    if (status) {
      fits_report_error(stderr,status);
      delete[] flux;
      return status;
    }

    fits_write_col(fp, TFLOAT, 1, 1, 1, nrows, flux,
                   &status);
    if (status) {
      fits_report_error(stderr,status);
      delete[] flux;
      return status;
    }
    
    fits_close_file(fp, &status);
    if (status) {
      fits_report_error(stderr,status);
      delete[] flux;
      return status;
    }

    delete[] flux;

  } catch ( const pofdExcept& ex ) {
    std::cout << "Error encountered" << std::endl;
    std::cout << ex << std::endl;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cout << "Bad allocation error: " << ba.what() << std::endl;
    return 16;
  } 

  return 0;
}

///////////////////////////////////

int makeCatDouble(int argc, char **argv) {

  unsigned int n0;
  std::string modelfile, outputfile; 
  unsigned long long int seed;
  bool verbose, have_user_seed;

  //Defaults
  verbose             = false;
  have_user_seed      = false;
  seed                = 0;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'S' :
      have_user_seed = true;
      seed = static_cast<unsigned long long int>( atoi(optarg) );
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc-2 ) {
    std::cout << "Some required arguments missing" << std::endl;
    std::cout << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  modelfile  = std::string(argv[optind]);
  n0         = atoi(argv[optind + 1]);
  outputfile = std::string(argv[optind + 2]);

  if (n0 == 0) {
    std::cout << "Invalid (non-positive) n0: " << n0 << std::endl;
    return 1;
  }

  try {

    //Set up random number generator
    ran rangen;
    if (!have_user_seed) {
      seed = static_cast<unsigned long long int>(time(NULL));
      seed += static_cast<unsigned long long int>(clock());
    }
    rangen.set_seed(seed);

    //Set up the model
    numberCountsDouble model(modelfile);
    if (!model.isValid())
      throw pofdExcept("pofd_coverage_makeCat", "makeCatDouble",
		       "Trying to realize model with invalid parameters", 1);

    //Get the fluxes
    float *flux1;
    float *flux2;
    flux1 = new float[n0];
    flux2 = new float[n0];
    std::pair<double, double> f12;
    for (unsigned int i = 0; i < n0; ++i) {
      f12 = model.genSource(rangen.doub(), rangen.gauss());
      flux1[i] = f12.first;
      flux2[i] = f12.second;
    }

    if (verbose) std::cout << "Writing catalog to " << outputfile 
			   << std::endl;

    //Do the write
    int status = 0;
    fitsfile *fp;
    fits_create_file(&fp, outputfile.c_str(), &status);
    if (status) {
      fits_report_error(stderr,status);
      return status;
    }

    int tfields = 2;
    long nrows = static_cast<long>(n0);
    char extname[] = "CATALOG";
    //These make g++ complain, but there isn't much I can do about it
    // because I can't change the call signature of fits_create_tbl
    // to make them const*
    char *ttype[] = {"FluxDensity1","FluxDensity2"};
    char *tform[] = {"1E","1E"};
    char *tunit[] = {"Jy","Jy"};
    fits_create_tbl(fp, BINARY_TBL, nrows, tfields, ttype, tform,
		    tunit, extname, &status);
    if (status) {
      fits_report_error(stderr,status);
      return status;
    }

    fits_write_col(fp, TFLOAT, 1, 1, 1, nrows, flux1,
                   &status);
    fits_write_col(fp, TFLOAT, 2, 1, 1, nrows, flux2,
                   &status);
    if (status) {
      fits_report_error(stderr,status);
      delete[] flux1; 
      delete[] flux2;
      return status;
    }
    
    fits_close_file(fp, &status);
    if (status) {
      fits_report_error(stderr,status);
      delete[] flux1; 
      delete[] flux2;
      return status;
    }

    delete[] flux1;
    delete[] flux2;

  } catch ( const pofdExcept& ex ) {
    std::cout << "Error encountered" << std::endl;
    std::cout << ex << std::endl;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cout << "Bad allocation error: " << ba.what() << std::endl;
    return 16;
  } 

  return 0;
}

///////////////////////////////////

int main( int argc, char** argv ) {
  bool twod;

  twod = false;

  //Only interested in a) displaying help and b) figuring out
  // if this is 1D or 2D c) displaying the version number
  int c;
  int option_index = 0;
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cout << "NAME" << std::endl;
      std::cout << "\tpofd_coverage_makeCat -- make simulated catalog for"
		<< " a spline type model." << std::endl;
      std::cout << std::endl;
      std::cout << "SYNOPSIS" << std::endl;
      std::cout << "\t  pofd_coverage_makeCat [options] modelfile n0 "
		<< "outputfile" << std::endl;
      std::cout << "DESCRIPTION" << std::endl;
      std::cout << "\tCreates a simulated catalog for a given model, and writes"
		<< " them" << std::endl;
      std::cout << "\tto outfile.  The number counts model in 1D is a spline"
		<< std::endl;
      std::cout << "\ta model specified by modelfile, and the number of"
		<< std::endl;
      std::cout << "\tsources to generate is n0." << std::endl;
      std::cout << std::endl;
      std::cout << "\tIn the 2D case the model is the 1D model in band 1 times"
		<< " a" << std::endl;
      std::cout << "\tLog-Normal distribution in flux2/flux1.  The mu and sigma"
		<< " Log-Normal" << std::endl;
      std::cout << "\tmodel parameters are stored as splines as a function of"
		<< " the" << std::endl;
      std::cout << "\tflux in the first band." << std::endl;
      std::cout << std::endl;
      std::cout << "\tmodelfile should be a text file.  In the 1D case it"
		<< " should" << std::endl;
      std::cout << "\tconsist of nknots lines of the form: " 
		<< std::endl << std::endl;
      std::cout << "\t\tflux_density n" << std::endl << std::endl;
      std::cout << "\twhere flux_density gives the positions of the knots"
		<< " in Jy" << std::endl;
      std::cout << "\tand n is the log10 differential number counts in"
		<< " deg^-2 Jy^-1" << std::endl;
      std::cout << "\tat the corresponding flux density.  Additional entries" 
		<< " on each"<< std::endl;
      std::cout << "\tline are ignored, and # denotes a comment line."
		<< std::endl;
      std::cout << std::endl;
      std::cout << "\tIn the 2D case the file should start with a line giving"
		<< " the" << std::endl;
      std::cout << "\tnumber of knots in the band 1 model, the number of"
		<< " knots in" << std::endl;
      std::cout << "\tthe sigma spline, and then the number in the mu spline."
		<< " This" << std::endl;
      std::cout << "\tshould be followed by nknots + nspline + nmu lines"
		<< std::endl;
      std::cout << "\tof the same form as the 1D model, with the first nknots"
		<< std::endl;
      std::cout << "\tspecifying the band 1 model as in the 1D case, and the"
		<< std::endl;
      std::cout << "\tfollowing lines giving the knot positions and values"
		<< " for" << std::endl;
      std::cout << "\tof the sigma and mu splines." << std::endl;
      std::cout << std::endl;
      std::cout << "OPTIONS" << std::endl;
      std::cout << "\t-d, --double" << std::endl;
      std::cout << "\t\tUse the 2D model." << std::endl;
      std::cout << "\t-h --help" << std::endl;
      std::cout << "\t\tPrint this message and exit." << std::endl;
      std::cout << "\t--S, --seed SEED" << std::endl;
      std::cout << "\t\tUse this seed for the random number generator." 
		<< std::endl;
      std::cout << "\t-v, --verbose" << std::endl;
      std::cout << "\t\tPrint informational messages while running"
		<< std::endl;
      std::cout << "\t-V, --version" << std::endl;
      std::cout << "\t\tOutput version number and exit" << std::endl;
      std::cout << "ONE-DIMENSIONAL OPTIONS" << std::endl;
      return 0;
      break;
    case 'd' :
      twod = true;
      break;
    case 'V' :
      std::cout << "pofd_coverage version number: " << pofd_coverage::version 
		<< std::endl;
      return 0;
      break;
    }

  if (twod)
    return makeCatDouble(argc, argv);
  else
    return makeCatSingle(argc, argv);
}



