#include<iostream>
#include<string>
#include<iomanip>
#include<vector>
#include<fstream>
#include<cstdlib>

#include<getopt.h>

#include<fitsio.h>
#include<fftw3.h>

#include "../include/global_settings.h"
#include "../include/beam.h"
#include "../include/doublebeam.h"
#include "../include/pofdExcept.h"
#include "../include/hipassFilter.h"

static struct option long_options[] = {
  {"double", no_argument, 0, 'd'},
  {"help", no_argument, 0, 'h'},
  {"histogram" , no_argument, 0, 'H'},
  {"inverse", no_argument, 0, 'i'},
  {"filterscale", required_argument, 0, 'F'},
  {"nfwhm", required_argument, 0, 'N'},
  {"nbins", required_argument, 0, '0'},
  {"oversamp", required_argument, 0, 'o'},
  {"verbose", no_argument, 0, 'v'},
  {"version", no_argument, 0, 'V'},
  {0,0,0,0}
};
char optstring[] = "dhHiF:N:0:o:vV";

///////////////////////////////

int getBeamSingle(int argc, char **argv) {

  unsigned int nbins, oversamp;
  double fwhm, nfwhm, pixsize, filterscale;
  bool verbose, histogram, inverse;
  std::string outputfile; // Output FITS file

  //Defaults
  nbins               = 120;
  nfwhm               = 3.5;
  verbose             = false;
  histogram           = false;
  oversamp            = 1;
  filterscale         = 0.0;
  inverse = false;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'F' :
      filterscale = atof(optarg);
      break;
    case 'H':
      histogram = true;
      break;
    case 'i':
      inverse = true;
      break;
    case 'N':
      nfwhm = atof(optarg);
      break;
    case '0':
      nbins = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'o':
      oversamp = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc - 2) {
    std::cout << "Some required arguments missing" << std::endl;
    std::cout << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  fwhm       = atof(argv[optind]);
  pixsize    = atof(argv[optind + 1]);
  outputfile = std::string(argv[optind + 2]);

  // Input checks
  if (histogram && (nbins == 0)) {
    std::cout << "Invalid (non-positive) number of beam histogram bins"
	      << std::endl;
    return 1;
  }
  if (fwhm <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM" << std::endl;
    return 1;
  }
  if (nfwhm <= 0) {
    std::cout << "Invalid (non-positive) number of beam FWHMs"
	      << std::endl;
    return 1;
  }
  if (pixsize >= fwhm / 2.0) {
    std::cout << "Insufficient (FWHM/2) beam sampling based on pixel size"
	      << std::endl;
    return 1;
  }
  if (oversamp % 2 == 0) {
    std::cout << "Invalid (non-odd) oversampling " << oversamp << std::endl;
    return 1;
  }
  if (filterscale < 0.0) {
    std::cout << "Invalid (negative) filter scale " << filterscale << std::endl;
    return 1;
  }

  try {
    beam bm(fwhm);
    
    if (verbose) {
      printf("   Beam fwhm:          %0.2f\n", bm.getFWHM());
      printf("   Beam area:          %0.3e\n", bm.getEffectiveArea());
      printf("   Pixel size:         %0.2f\n", pixsize);
      if (filterscale > 0.0)
	printf("   filter scale:       %0.4f\n", filterscale);
      if (oversamp != 1)
	printf("oversamp:              %u\n", oversamp);
      if (inverse) printf("Returning inverse beam\n");
      if (histogram) printf("Returning histogrammed beam\n");
    }

    if (histogram) {
      // Get histogrammed beam
      beamHist bmhist(nbins, filterscale);
      bmhist.fill(bm, nfwhm, pixsize, inverse, oversamp);
      // Write
      bmhist.writeToFits(outputfile);
    } else {
      hipassFilter *filt = NULL;
      if (filterscale > 0) filt = new hipassFilter(filterscale);
      bm.writeToFits(outputfile, pixsize, nfwhm, oversamp, filt, inverse);
      if (filt != NULL) delete filt; 
    }    
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

///////////////////////////////////////

int getBeamDouble(int argc, char **argv) {

  bool verbose, histogram, inverse;
  unsigned int nbins, oversamp;
  double fwhm1, fwhm2, nfwhm, pixsize, filterscale;
  std::string outputfile; //Ouput pofd option

  //Defaults
  nbins               = 150;
  nfwhm               = 3.5;
  verbose             = false;
  histogram           = false;
  filterscale         = 0.0;
  oversamp            = 1;
  inverse             = false;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'F' :
      filterscale = atof(optarg);
      break;
    case 'H':
      histogram = true;
      break;
    case 'i':
      inverse = true;
      break;
    case 'N':
      nfwhm = atof(optarg);
      break;
    case '0':
      nbins = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'o':
      oversamp = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc - 3) {
    std::cout << "Some required arguments missing" << std::endl;
    std::cout << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  fwhm1      = atof(argv[optind]);
  fwhm2      = atof(argv[optind + 1]);
  pixsize    = atof(argv[optind + 2]);
  outputfile = std::string(argv[optind + 3]);

  //Input tests
  if (histogram && (nbins == 0)) {
    std::cout << "Invalid (non-positive) number of beam histogram bins"
	      << std::endl;
    return 1;
  }
  if (fwhm1 <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM (band 1)" << std::endl;
    return 1;
  }
  if (fwhm2 <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM (band 2)" << std::endl;
    return 1;
  }
  if (nfwhm <= 0) {
    std::cout << "Invalid (non-positive) number of beam FWHMs"
	      << std::endl;
    return 1;
  }
  if (pixsize >= fwhm1 / 2.0 || pixsize >= fwhm2 / 2.0) {
    std::cout << "Insufficient (FWHM/2) beam sampling based on pixel size"
	      << std::endl;
    return 1;
  }

  try {
    doublebeam bm(fwhm1, fwhm2);

    if (verbose) {
      std::pair<double, double> dpr;
      dpr = bm.getFWHM();
      printf("   Beam fwhm1:         %0.2f\n", dpr.first);
      printf("   Beam fwhm2:         %0.2f\n", dpr.second);
      dpr = bm.getEffectiveArea();
      printf("   Beam area1:         %0.3e\n", dpr.first);
      printf("   Beam area2:         %0.3e\n", dpr.second);
      printf("   Pixel size:         %0.2f\n", pixsize);
      if (filterscale > 0.0)
	printf("   filter scale:       %0.4f\n", filterscale);
      if (oversamp != 1)
	printf("oversamp:              %u\n", oversamp);
      if (inverse) printf("Returning inverse beam\n");
      if (histogram) printf("Returning histogrammed beam\n");
    }

    if (histogram) {
      // Get histogrammed beam
      doublebeamHist bmhist(nbins, filterscale);
      bmhist.fill(bm, nfwhm, pixsize, inverse, oversamp);
      // Write
      bmhist.writeToFits(outputfile);
    } else {
      hipassFilter *filt = NULL;
      if (filterscale > 0) filt = new hipassFilter(filterscale);
      bm.writeToFits(outputfile, pixsize, nfwhm, oversamp, filt, inverse);
      if (filt != NULL) delete filt; 
    }    
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


////////////////////////////////////
int main( int argc, char** argv ) {
  bool twod;

  twod = false;

  //Only interested in a) displaying help and b) figuring out
  // if this is 1D or 2D c) displaying the version number
  int c;
  int option_index = 0;
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'h' :
      std::cout << "NAME" << std::endl;
      std::cout << "\tpofd_coverage_getBeam -- get the beam for the"
		<< std::endl;
      std::cout << "\t pofd_coverage model." << std::endl;
      std::cout << std::endl;
      std::cout << "SYNOPSIS" << std::endl;
      std::cout << "\t One-dimensional case:" << std::endl;
      std::cout << "\t  pofd_coverage_getBeam [options] fwhm pixsize"
		<< " outfile" << std::endl; 
      std::cout << std::endl;
      std::cout << "\t Two-dimensional case:" << std::endl;
      std::cout << "\t  pofd_coverage_getPD [options] -d fwhm1 fwhm2 pixsize"
		<< std::endl;
      std::cout << "\t    outfile" << std::endl; 
      std::cout << std::endl;
      std::cout << "DESCRIPTION" << std::endl;
      std::cout << "\tWrites the beam out to a FITS file." << std::endl;
      std::cout << std::endl;
      std::cout << "OPTIONS" << std::endl;
      std::cout << "\t-d, --twod" << std::endl;
      std::cout << "\t\tIf set, the two-dimensional model is used."
		<< std::endl;
      std::cout << "\t-F, --filtscale VALUE" << std::endl;
      std::cout << "\t\tRadius of high-pass filter in arcseconds. If zero,"
		<< std::endl;
      std::cout << "\t\tno filtering is applied (def: 0)." << std::endl;
      std::cout << "\t-h --help" << std::endl;
      std::cout << "\t\tPrint this message and exit." << std::endl;
      std::cout << "\t-H --histogram" << std::endl;
      std::cout << "\t\tWrite out the histogrammed beam." << std::endl;
      std::cerr << "\t-i, --inverse" << std::endl;
      std::cout << "\t\tReturn the inverse beam rather than the beam."
		<< std::endl;
      std::cout << "\t-N, --nfwhm value" << std::endl;
      std::cout << "\t\tNumber of beam FWHM out to go when computing beam."
		<< "(def: 3.5)" << std::endl;
      std::cout << "\t--nbins value" << std::endl;
      std::cout << "\t\tNumber of bins to use in histogrammed beam. (def: 120)"
		<< std::endl;
      std::cout << "\t-o, --oversample VALUE" << std::endl;
      std::cout << "\t\tAmount to oversample the beam; must be odd integer."
		<< " (def: 1)" << std::endl;
      std::cout << "\t-v, --verbose" << std::endl;
      std::cout << "\t\tPrint informational messages while running"
		<< std::endl;
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

  if (!twod)
    return getBeamSingle(argc, argv);
  else
    return getBeamDouble(argc, argv);
}
