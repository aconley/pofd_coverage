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
#include "../include/fourierFilter.h"

static struct option long_options[] = {
  {"double", no_argument, 0, 'd'},
  {"help", no_argument, 0, 'h'},
  {"histogram" , no_argument, 0, 'H'},
  {"inverse", no_argument, 0, 'i'},
  {"filterscale", required_argument, 0, 'F'},
  {"matched", no_argument, 0, 'm'},
  {"nfwhm", required_argument, 0, 'N'},
  {"nkeep", required_argument, 0, '1'},
  {"nbins", required_argument, 0, '0'},
  {"oversamp", required_argument, 0, 'o'},
  {"sigc", required_argument, 0, '3'},
  {"sigi", required_argument, 0, '2'},
  {"sigi1", required_argument, 0, '4'},
  {"sigi2", required_argument, 0, '5'},
  {"verbose", no_argument, 0, 'v'},
  {"version", no_argument, 0, 'V'},
  {0, 0, 0, 0}
};
char optstring[] = "dhHiF:mN:1:0:o:2:3:4:5:vV";

///////////////////////////////

int getBeamSingle(int argc, char **argv) {

  unsigned int nbins, oversamp;
  bool verbose, histogram, inverse, matched;
  double fwhm, nfwhm, pixsize, filterscale, nkeep, sigi, sigc;
  std::string outputfile; // Output FITS file

  //Defaults
  nbins               = 120;
  nfwhm               = 3.5;
  nkeep               = 0.0;
  verbose             = false;
  histogram           = false;
  oversamp            = 1;
  filterscale         = 0.0;
  matched             = false;
  sigi                = 0.002;
  sigc                = 0.006;
  inverse             = false;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'F':
      filterscale = atof(optarg);
      break;
    case 'H':
      histogram = true;
      break;
    case 'i':
      inverse = true;
      break;
    case 'm':
      matched = true;
      break;
    case 'N':
      nfwhm = atof(optarg);
      break;
    case '1':
      nkeep = atof(optarg);
      break;
    case '0':
      nbins = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'o':
      oversamp = static_cast<unsigned int>(atoi(optarg));
      break;
    case '2':
      sigi = atof(optarg);
      break;
    case '3':
      sigc = atof(optarg);
      break;
    case 'v':
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
  if (nkeep < 0) {
    std::cout << "Invalid (negative) nkeep: " << nkeep << std::endl;
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
  if (matched) {
    if (sigi <= 0.0) {
      std::cout << "Invalid (non-positive) sigi for matched filter"
		<< std::endl;
      return 1;
    }
    if (sigc <= 0.0) {
      std::cout << "Invalid (non-positive) sigc for matched filter"
		<< std::endl;
      return 1;
    }
  }

  try {
    beam bm(fwhm);
    
    if (verbose) {
      printf("   Beam fwhm:          %0.1f\n", bm.getFWHM());
      printf("   Beam area:          %0.3e\n", bm.getEffectiveArea());
      printf("  Pixel size:          %0.2f\n", pixsize);
      if (filterscale > 0.0)
	printf("filter scale:          %0.2f\n", filterscale);
      if (matched) {
	printf("matched fwhm:          %0.1f\n", bm.getFWHM());
	printf("matched sigi:          %0.4f\n", sigi);
	printf("matched sigc:          %0.4f\n", sigc);
      }
      if (oversamp != 1)
	printf("    oversamp:          %u\n", oversamp);
      if (inverse) printf("Returning inverse beam\n");
      if (histogram) printf("Returning histogrammed beam\n");
    }

    fourierFilter *filt = NULL;
    if (filterscale > 0) {
      if (matched)
	filt = new fourierFilter(pixsize, fwhm, sigi, sigc,
				 filterscale, 0.1, true);
      else
	filt = new fourierFilter(pixsize, filterscale, 0.1, true);
    } else if (matched)
	filt = new fourierFilter(pixsize, fwhm, sigi, sigc, true);

    if (histogram) {
      // Get histogrammed beam
      beamHist bmhist(nbins);
      bmhist.fill(bm, nfwhm, pixsize, inverse, oversamp, filt, nkeep);

      // Write
      bmhist.writeToFits(outputfile);
    } else 
      bm.writeToFits(outputfile, pixsize, nfwhm, oversamp, filt, inverse);

    if (filt != NULL) delete filt;
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

  bool verbose, histogram, inverse, matched;
  unsigned int nbins, oversamp;
  double fwhm1, fwhm2, nfwhm, pixsize, nkeep;
  double filterscale, sigi1, sigi2, sigc; // Filtering params
  std::string outputfile; //Ouput pofd option

  //Defaults
  nbins               = 150;
  nfwhm               = 3.5;
  verbose             = false;
  histogram           = false;
  filterscale         = 0.0;
  matched             = false;
  sigc                = 0.006;
  sigi1               = 0.002;
  sigi2               = 0.002;
  oversamp            = 1;
  inverse             = false;
  nkeep               = 0;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'F':
      filterscale = atof(optarg);
      break;
    case 'H':
      histogram = true;
      break;
    case 'i':
      inverse = true;
      break;
    case 'm':
      matched = true;
      break;
    case 'N':
      nfwhm = atof(optarg);
      break;
    case '1':
      nkeep = atof(optarg);
      break;
    case '0':
      nbins = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'o':
      oversamp = static_cast<unsigned int>(atoi(optarg));
      break;
    case '3':
      sigc = atof(optarg);
      break;
    case '4':
      sigi1 = atof(optarg);
      break;
    case '5':
      sigi2 = atof(optarg);
      break;
    case 'v':
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
  if (nkeep < 0) {
    std::cout << "Invalid (negative) nkeep: " << nkeep << std::endl;
    return 1;
  }
  if (pixsize >= fwhm1 / 2.0 || pixsize >= fwhm2 / 2.0) {
    std::cout << "Insufficient (FWHM/2) beam sampling based on pixel size"
	      << std::endl;
    return 1;
  }
  if (matched) {
    if (sigi1 <= 0.0) {
      std::cout << "Invalid (non-positive) sigi1 for matched filter"
		<< std::endl;
      return 1;
    }
    if (sigi2 <= 0.0) {
      std::cout << "Invalid (non-positive) sigi2 for matched filter"
		<< std::endl;
      return 1;
    }
    if (sigc <= 0.0) {
      std::cout << "Invalid (non-positive) sigc for matched filter"
		<< std::endl;
      return 1;
    }
  }

  try {
    doublebeam bm(fwhm1, fwhm2);

    if (verbose) {
      std::pair<double, double> dpr;
      dpr = bm.getFWHM();
      printf("   Beam fwhm1:         %0.1f\n", dpr.first);
      printf("   Beam fwhm2:         %0.1f\n", dpr.second);
      dpr = bm.getEffectiveArea();
      printf("   Beam area1:         %0.3e\n", dpr.first);
      printf("   Beam area2:         %0.3e\n", dpr.second);
      printf("   Pixel size:         %0.2f\n", pixsize);
      if (filterscale > 0.0)
	printf(" filter scale:          %0.1f\n", filterscale);
      if (matched) {
	printf(" matched fwhm1:         %0.1f\n", dpr.first);
	printf(" matched fwhm2:         %0.1f\n", dpr.second);
	printf(" matched sigi1:         %0.4f\n", sigi1);
	printf(" matched sigi2:         %0.4f\n", sigi2);
	printf(" matched sigc:          %0.4f\n", sigc);
      }
      if (oversamp != 1)
	printf("    oversamp:          %u\n", oversamp);
      if (inverse) printf("Returning inverse beam\n");
      if (histogram) printf("Returning histogrammed beam\n");
    }

    fourierFilter *filt1 = NULL, *filt2 = NULL;
    if (filterscale > 0) {
      if (matched) {
	filt1 = new fourierFilter(pixsize, fwhm1, sigi1, sigc,
				  filterscale, 0.1, true);
	filt2 = new fourierFilter(pixsize, fwhm2, sigi2, sigc,
				  filterscale, 0.1, true);
      } else
	filt1 = new fourierFilter(pixsize, filterscale, 0.1, true);
    } else if (matched) {
      filt1 = new fourierFilter(pixsize, fwhm1, sigi1, sigc, true);
      filt2 = new fourierFilter(pixsize, fwhm1, sigi2, sigc, true);
    }
    if (histogram) {
      // Get histogrammed beam
      doublebeamHist bmhist(nbins);
      bmhist.fill(bm, nfwhm, pixsize, inverse, oversamp, filt1, filt2, nkeep);
      // Write
      bmhist.writeToFits(outputfile);
    } else
      bm.writeToFits(outputfile, pixsize, nfwhm, oversamp, filt1, filt2, 
		     inverse);

    if (filt1 != NULL) delete filt1; 
    if (filt2 != NULL) delete filt2; 

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
    case 'h':
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
      std::cout << "\t-m, --matched" << std::endl;
      std::cout << "\t\tApply matched filtering to the beam, with a FWHM"
		<< " matching the" << std::endl;
      std::cout << "\t\tbeam (the band 1 beam in the 2d case) and the "
		<< " instrument" << std::endl;
      std::cout << "\t\tand confusion noise controlled by --sigi and --sigc,"
		<< " or" << std::endl;
      std::cout << "\t\t--sigi1, --sigi2 and --sigc in 2D." << std::endl;
      std::cout << "\t\tOff by default." << std::endl;
      std::cout << "\t--nbins value" << std::endl;
      std::cout << "\t\tNumber of bins to use in histogrammed beam. (def: 120)"
		<< std::endl;
      std::cout << "\t-N, --nfwhm VALUE" << std::endl;
      std::cout << "\t\tNumber of beam FWHM out to go when computing beam."
		<< "(def: 3.5)" << std::endl;
      std::cout << "\t--nkeep VALUE" << std::endl;
      std::cout << "\t\tNumber of FWHM to keep after histogramming.  Only"
		<< " applies" << std::endl;
      std::cout << "\t\tif the beam is histogrammed.  The default is to keep"
		<< std::endl;
      std::cout << "\t\tall of the beam specified by --nfwhm." << std::endl;
      std::cout << "\t-o, --oversample VALUE" << std::endl;
      std::cout << "\t\tAmount to oversample the beam; must be odd integer."
		<< " (def: 1)" << std::endl;
      std::cout << "\t--sigc VALUE" << std::endl;
      std::cout << "\t\tConfusion noise for matched filtering, in Jy. (Def:"
		<< " 0.006)" << std::endl;
      std::cout << "\t-v, --verbose" << std::endl;
      std::cout << "\t\tPrint informational messages while running"
		<< std::endl;
      std::cout << "ONE-DIMENSIONAL OPTIONS" << std::endl;
      std::cout << "\t--sigi VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, in Jy. (Def:"
		<< " 0.002)" << std::endl;
      std::cout << "TWO-DIMENSIONAL OPTIONS" << std::endl;
      std::cout << "\t--sigi1 VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, band 1, in "
		<< std::endl;
      std::cout << "\t\tJy. (Def: 0.002)" << std::endl;
      std::cout << "\t--sigi2 VALUE" << std::endl;
      std::cout << "\t\tInstrument noise for matched filtering, band 2, in "
		<< std::endl;
      std::cout << "\t\tJy. (Def: 0.002)." << std::endl;
      std::cout << std::endl;
      return 0;
      break;
    case 'd':
      twod = true;
      break;
    case 'V':
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
