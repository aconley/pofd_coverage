#include<iostream>
#include<string>
#include<iomanip>
#include<vector>
#include<fstream>
#include<cstdlib>

#include<getopt.h>

#include<global_settings.h>
#include<utility.h>
#include<PDFactory.h>
#include<beam.h>
#include<numberCounts.h>
#include<paramSet.h>
#include<pofdExcept.h>

static struct option long_options[] = {
  {"double", no_argument, 0, 'd'},
  {"edgeinterp",required_argument,0,'e'},
  {"help", no_argument, 0, 'h'},
  {"fits", no_argument, 0, 'f'},
  {"log", no_argument, 0, 'l'},
  {"nflux", required_argument, 0, 'n'},
  {"nfwhm", required_argument, 0, 'N'},
  {"nbins", required_argument, 0, '0'},
  {"pixsize", required_argument, 0, 'p'},
  {"rfile", required_argument, 0, 'r'},
  {"sigma", required_argument, 0, 's'},
  {"verbose", no_argument, 0, 'v'},
  {"version", no_argument, 0, 'V'},
  {"wisdom", required_argument, 0, 'w'},
  {0,0,0,0}
};
char optstring[] = "dhe:fln:N:0:p:r:s:vVw:";

///////////////////////////////

int getPDSingle(int argc, char **argv) {

  std::string modelfile; //Knot parameters
  double n0; //Number of sources
  double maxflux, fwhm, nfwhm, pixsize; //Calculation params (req)
  double sigma; //Instrument noise
  std::string outputfile; //Ouput pofd option
  unsigned int nflux, nbins;
  bool has_wisdom, verbose, return_log, write_fits, has_user_pixsize, write_r;
  std::string wisdom_file, r_file;

  //Defaults
  sigma               = 2e-3;
  has_wisdom          = false;
  nflux               = 131072;
  nbins               = 80;
  nfwhm               = 3.5;
  verbose             = false;
  return_log          = false;
  write_fits          = false;
  has_user_pixsize    = false;
  write_r             = false;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'f' :
      write_fits = true;
      break;
    case 'l' :
      return_log = true;
      break;
    case 'n' :
      nflux = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'N':
      nfwhm = atof(optarg);
      break;
    case '0':
      nbins = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'p':
      has_user_pixsize = true;
      pixsize = atof(optarg);
      break;
    case 'r':
      write_r = true;
      r_file = std::string(optarg);
    case 's' :
      sigma = atof(optarg);
      break;
    case 'v' :
      verbose = true;
      break;
      return 0;
      break;
    case 'w' :
      has_wisdom = true;
      wisdom_file = std::string(optarg);
      break;
    }

  if (optind >= argc - 4) {
    std::cerr << "Some required arguments missing" << std::endl;
    std::cerr << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  modelfile  = std::string(argv[optind]);
  n0         = atof(argv[optind + 1]);
  fwhm       = atof(argv[optind + 2]);
  maxflux    = atof(argv[optind + 3]);
  outputfile = std::string(argv[optind + 4]);

  if (!has_user_pixsize) pixsize = fwhm / 3.0;

  //Input tests
  if (nflux == 0) {
    std::cerr << "Error -- number of fluxes requested is zero."
	      << std::endl;
    return 1;
  }
  if (nflux & (nflux-1)) {
    std::cerr << "nflux must be power of 2" << std::endl;
    std::cerr << " yours is: " << nflux << std::endl;
    return 1;
  }
  if (sigma < 0.0) {
    std::cerr << "Invalid noise level: must be >= 0.0" << std::endl;
    return 1;
  }
  if (n0 <= 0.0) {
    std::cerr << "Invalid (non-positive) number of sources per area"
	      << std::endl;
    return 1;
  }
  if (nbins == 0) {
    std::cerr << "Invalid (non-positive) number of beam histogram bins"
	      << std::endl;
    return 1;
  }
  if (fwhm <= 0.0) {
    std::cerr << "Invalid (non-positive) FWHM" << std::endl;
    return 1;
  }
  if (nfwhm <= 0) {
    std::cerr << "Invalid (non-positive) number of beam FWHMs"
	      << std::endl;
    return 1;
  }
  if (pixsize >= fwhm/2.0) {
    std::cerr << "Insufficient (FWHM/2) beam sampling based on pixel size"
	      << std::endl;
    return 1;
  }

  try {
    numberCounts model(modelfile);
    beam bm(fwhm);
    PDFactory pfactory;
    PD pd;

    if (verbose) pfactory.setVerbose();

    bool success;
    if (has_wisdom) {
      std::cout << "Reading in wisdom file: " << wisdom_file 
		<< std::endl;
      success = pfactory.addWisdom(wisdom_file);
      if (!success) {
	std::cerr << "Error reading wisdom file: " << wisdom_file << std::endl;
	return 4;
      }
    }
    
    if (verbose) {
      printf("   Beam fwhm:          %0.2f\n", bm.getFWHM());
      printf("   Beam area:          %0.3e\n", bm.getEffectiveArea());
      printf("   Mean flux per area: %0.2f\n",
	     model.getMeanFluxPerArea());
      printf("   Base N0:            %0.4e\n", model.getBaseN0());
      printf("   N0:                 %0.4e\n", n0);
      printf("sigma:                 %0.4f\n", sigma);
      if (return_log) 
	printf("  Returning log( P(D) ) rather than P(D)\n");
    }

    //Get P(D)
    if (verbose) std::cout << "Getting P(D) with transform length: " 
			   << nflux << " and max flux: " 
			   << maxflux << std::endl;
    pfactory.initPD(nflux, sigma, maxflux, n0, model, bm,
		    pixsize, nfwhm, nbins);
    pfactory.getPD(n0, pd, return_log, true);
    
   if (write_r) {
      if (verbose) std::cout << "Writing R to " << r_file << std::endl;
      pfactory.writeRToFile(r_file);
    }
    
    //Write it
    if (verbose) std::cout << "Writing P(D) to " << outputfile 
			   << std::endl;
    if (verbose && write_fits)
      std::cout << " Writing as FITS file" << std::endl;
    if (write_fits) {
      pd.writeToFits(outputfile);
    } else {
      std::ofstream ofs(outputfile.c_str());
      if (!ofs) {
	std::cerr << "Error opening output file: " << outputfile
		  << std::endl;
	return 64;
      }
      ofs << pd;
      ofs.close();
    }
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

///////////////////////////////////////


////////////////////////////////////
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
      std::cerr << "\tpofd_coverage_getPD -- get the P(D) for a broken power"
		<< " law" << std::endl;
      std::cerr << "\t type model with a Gaussian beam." << std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\t  pofd_coverage_getPD [options] modelfile n0 fwhm maxflux"
		<< " outfile" << std::endl; 
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEvaluates P(D) for the specified model and writes it to" 
		<< std::endl;
      std::cerr << "\toutfile.  The number counts model is a broken power" 
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
      std::cerr << "\tand n is the log10 differential number counts in"
		<< " deg^-2 Jy^-1" << std::endl;
      std::cerr << "\tat the corresponding flux density.  Additional entries" 
		<< " on each"<< std::endl;
      std::cerr << "\tline are ignored, and # denotes a comment line."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tfwhm is the beam FWHM in arcsec.  The beam is assumed "
		<< "Gaussian. " << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tmaxflux is the maximum flux density generated.  In"
		<< " general" << std::endl;
      std::cerr << "\tthe maxflux values will not quite be realized."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-d, --twod" << std::endl;
      std::cerr << "\t\tIf set, the two-dimensional model is used."
		<< std::endl;
      std::cerr << "\t-h --help" << std::endl;
      std::cerr << "\t\tPrint this message and exit." << std::endl;
      std::cerr << "\t-f, --fits" << std::endl;
      std::cerr << "\t\tWrite output as a fits file rather than text."
		<< std::endl;
      std::cerr << "\t-l, --log" << std::endl;
      std::cerr << "\t\tReturn the log P(D) rather than the P(D)."
		<< std::endl;
      std::cerr << "\t-n, --nflux value" << std::endl;
      std::cerr << "\t\tThe number of requested fluxes along each dimension."
		<< std::endl;
      std::cerr << "\t\tAlso sets the transform size used. (def: 131072)."
		<< std::endl;
      std::cerr << "\t-N, --nfwhm value" << std::endl;
      std::cerr << "\t\tNumber of beam FWHM out to go when computing beam."
		<< "(def: 3.5)" << std::endl;
      std::cerr << "\t--nbins value" << std::endl;
      std::cerr << "\t\tNumber of bins to use in histogrammed beam. (def: 80)"
		<< std::endl;
      std::cerr << "\t-p, --pixsize value" << std::endl;
      std::cerr << "\t\tPixel size in arcsec. (def: FWHM/3.0)" << std::endl;
      std::cerr << "\t-r, --rfile FILENAME" << std::endl;
      std::cerr << "\t\tWrite the R used to this file as text." << std::endl;
      std::cerr << "\t-s, --sigma VALUE" << std::endl;
      std::cerr << "\t\tThe assumed per-pixel noise (def: 0.002)" << std::endl;
      std::cerr << "\t-v, --verbose" << std::endl;
      std::cerr << "\t\tPrint informational messages while running"
		<< std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput version number and exit" << std::endl;
      std::cerr << "\t-w, --wisdom wisdomfile" << std::endl;
      std::cerr << "\t\tName of wisdom file (prepared with fftw-wisdom)." 
		<< std::endl;
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
    return getPDSingle(argc,argv);
  else {
    std::cerr << "2D model not supported" << std::endl;
    return 1;
  }
}
