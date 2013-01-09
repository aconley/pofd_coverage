#include<iostream>
#include<string>
#include<iomanip>
#include<vector>
#include<fstream>
#include<cstdlib>

#include<getopt.h>

#include<simImage.h>
#include<pofdExcept.h>
#include<global_settings.h>

//Set up global option index that can be used for both single and double case
static struct option long_options[] = {
  {"double",no_argument,0,'d'},
  {"extra_smooth",required_argument,0,'e'},
  {"extra_smooth1",required_argument,0,'1'},
  {"extra_smooth2",required_argument,0,'2'},
  {"help", no_argument, 0, 'h'},
  {"oversample", required_argument, 0, 'o'},
  {"seed", required_argument, 0,'S'},
  {"sigma",required_argument,0,'s'},
  {"sigma1",required_argument,0,'3'},
  {"sigma2",required_argument,0,'4'},
  {"verbose",no_argument,0,'v'},
  {"version",no_argument,0,'V'},
  {0,0,0,0}
};
char optstring[] = "de:1:2:ho:S:s:3:4:vV";


int makeSimSingle(int argc, char **argv) {

  unsigned int n1, n2;
  double n0, pixsize, sigma, fwhm;
  double extra_smooth; //Additional smoothing
  std::string modelfile, outputfile; 
  unsigned long long int user_seed;
  bool verbose, do_extra_smooth, have_user_seed;
  unsigned int oversample;

  //Defaults
  do_extra_smooth     = false;
  extra_smooth        = 0.0;
  sigma               = 0.0;
  verbose             = false;
  user_seed           = 0;
  have_user_seed      = false;
  oversample          = 1;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'e' :
      do_extra_smooth = true;
      extra_smooth = atof(optarg);
      break;
    case 'o':
      oversample = atoi(optarg);
      break;
    case 'S' :
      have_user_seed = true;
      user_seed = static_cast<unsigned long long int>( atoi(optarg) );
      break;
    case 's' :
      sigma = atof(optarg);
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc-6 ) {
    std::cerr << "Some required arguments missing" << std::endl;
    std::cerr << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  modelfile  = std::string( argv[optind] );
  n0         = atof(argv[optind+1]);
  fwhm       = atof(argv[optind+2]);
  pixsize    = atof(argv[optind+3]);
  n1         = atoi(argv[optind+4]);
  n2         = atoi(argv[optind+5]);
  outputfile = std::string(argv[optind+6]);

  if (n0 <= 0.0) {
    std::cerr << "Invalid (non-positive) n0: " << n0 << std::endl;
    return 1;
  }
  if (std::isnan(n0) || std::isinf(n0)) {
    std::cerr << "Invalid (non-finite) n0: " << n0 << std::endl;
    return 1;
  }
  if (sigma < 0.0) {
    std::cerr << "Invalid noise level: must be >= 0.0" << std::endl;
    return 1;
  }
  if (fwhm < 0.0) {
    std::cerr << "Invalid (non-positive) FWHM" << std::endl;
    return 1;
  }
  if (do_extra_smooth && extra_smooth < 0.0) {
    std::cerr << "Invalid (non-positive) extra smoothing FWHM" << std::endl;
    return 1;
  }
  if (oversample == 0) {
    std::cerr << "Invalid (non-positive) oversampling" << std::endl;
    return 1;
  }
  

  try {
    numberCounts model(modelfile);
    simImage dim(n1, n2, pixsize, fwhm, sigma, extra_smooth,
		 oversample);
    if (have_user_seed) dim.setSeed( user_seed );
    dim.realize(model, n0, do_extra_smooth, true, false); //Do mean subtract

    //Write it
    if (verbose) std::cout << "Writing simulated image to " << outputfile 
			   << std::endl;
    int status = dim.writeToFits(outputfile);
    if (status != 0) return status;
  } catch ( const pofdExcept& ex ) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
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
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tpofd_coverage_makeSim -- make simulated images for"
		<< " a broken" << std::endl;
      std::cerr << "\t power law type model with Gaussian beams."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\t  pofd_coverage_makeSim [options] modelfile n0 fwhm "
		<< "pixsize" << std::endl;
      std::cerr << "\t    n1 n2 outputfile" << std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tCreates a simulated image for a given model, and writes"
		<< " them" << std::endl;
      std::cerr << "\tto outfile.  The number counts model is a broken power" 
		<< std::endl;
      std::cerr << "\ta law model specified by modelfile, and by the number of"
		<< std::endl;
      std::cerr << "\tsources per unit area n0." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tmodelfile should be a text file that consists of"
		<< " nknots" << std::endl;
      std::cerr << "\tlines of the form:" << std::endl << std::endl;
      std::cerr << "\t\tflux_density n" << std::endl << std::endl;
      std::cerr << "\twhere flux_density gives the positions of the knots"
		<< " in Jy" << std::endl;
      std::cerr << "\tand n is the differential number counts in deg^-2 Jy^-1"
		<< " at" << std::endl;
      std::cerr << "\tthe corresponding flux density.  Additional entries on" 
		<< " each"<< std::endl;
      std::cerr << "\tline are ignored, and # denotes a comment line."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tfwhm is the beam FWHM in arcsec.  The beam is assumed "
		<< "Gaussian. " << std::endl;
      std::cerr << "\tPixsize gives the pixel size (in arcsec), while n1 and "
		<< "n2 are the" << std::endl;
      std::cerr << "\tnumber of pixels along each dimension." << std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-h --help" << std::endl;
      std::cerr << "\t\tPrint this message and exit." << std::endl;
      std::cerr << "\t-e, --extra_smooth FWHM" << std::endl;
      std::cerr << "\t\tApply additional smoothing with a gaussian of this"
		<< " FWHM" << std::endl;
      std::cerr << "\t\t(in arcseconds); this is applied after noise is added."
		<< std::endl;
      std::cerr << "\t-o, --oversample VALUE" << std::endl;
      std::cerr << "\t\tAmount of oversampling to use (integral) when " 
		<< "generating" << std::endl;
      std::cerr << "\t\timage.  The data is then down-binned to the specified"
		<< "size." << std::endl;
      std::cerr << "\t\tThe default is to apply no oversampling." << std::endl;
      std::cerr << "\t--sigma NOISE" << std::endl;
      std::cerr << "\t\tThe assumed per-pixel noise (def: 0)." << std::endl;
      std::cerr << "\t--S, --seed SEED" << std::endl;
      std::cerr << "\t\tUse this seed for the random number generator." 
		<< std::endl;
      std::cerr << "\t-v, --verbose" << std::endl;
      std::cerr << "\t\tPrint informational messages while running"
		<< std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput version number and exit" << std::endl;
      std::cerr << "ONE-DIMENSIONAL OPTIONS" << std::endl;
      return 0;
      break;
    case 'd' :
      twod = true;
      break;
    case 'V' :
      std::cerr << "pofd_coverage version number: " << pofd_coverage::version 
		<< std::endl;
      return 0;
      break;
    }

  if (!twod)
    return makeSimSingle(argc,argv);
  else {
    std::cerr << "2D broken power law model not supported" << std::endl;
    return 1;
  }

}



