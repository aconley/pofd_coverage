#include<iostream>
#include<string>
#include<iomanip>
#include<vector>
#include<fstream>
#include<cstdlib>

#include<getopt.h>

#include "../include/global_settings.h"
#include "../include/utility.h"
#include "../include/PDFactory.h"
#include "../include/beam.h"
#include "../include/numberCounts.h"
#include "../include/PDFactoryDouble.h"
#include "../include/doublebeam.h"
#include "../include/numberCountsDouble.h"
#include "../include/pofdExcept.h"

static struct option long_options[] = {
  {"double", no_argument, 0, 'd'},
  {"edgeinterp",required_argument,0,'e'},
  {"help", no_argument, 0, 'h'},
  {"fits", no_argument, 0, 'f'},
  {"filterscale", required_argument, 0, 'F'},
  {"log", no_argument, 0, 'l'},
  {"nflux", required_argument, 0, 'n'},
  {"nfwhm", required_argument, 0, 'N'},
  {"nbins", required_argument, 0, '0'},
  {"oversamp", required_argument, 0, 'o'},
  {"pixsize", required_argument, 0, 'p'},
  {"rfile", required_argument, 0, 'r'},
  {"sigma", required_argument, 0, 's'},
  {"sigma1",required_argument,0,'3'},
  {"sigma2",required_argument,0,'4'},
  {"verbose", no_argument, 0, 'v'},
  {"version", no_argument, 0, 'V'},
  {"wisdom", required_argument, 0, 'w'},
  {0,0,0,0}
};
char optstring[] = "dhe:fF:ln:N:0:o:p:r:s:3:4:vVw:";

///////////////////////////////

int getPDSingle(int argc, char **argv) {

  std::string modelfile; //Knot parameters
  double n0; //Number of sources
  double maxflux, fwhm, nfwhm, pixsize, filterscale; //Calculation params (req)
  double sigma; //Instrument noise
  std::string outputfile; //Ouput pofd option
  unsigned int nflux, nbins;
  bool has_wisdom, verbose, return_log, write_fits, has_user_pixsize, write_r;
  std::string wisdom_file, r_file;
  unsigned int oversample;

  //Defaults
  sigma               = 2e-3;
  has_wisdom          = false;
  nflux               = 131072;
  nbins               = 120;
  nfwhm               = 4.5;
  verbose             = false;
  return_log          = false;
  write_fits          = false;
  has_user_pixsize    = false;
  pixsize             = 3.0;
  write_r             = false;
  oversample          = 1;
  filterscale         = 0.0;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'f' :
      write_fits = true;
      break;
    case 'F' :
      filterscale = atof(optarg);
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
    case 'o':
      oversample = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'p':
      has_user_pixsize = true;
      pixsize = atof(optarg);
      break;
    case 'r':
      write_r = true;
      r_file = std::string(optarg);
      break;
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
    std::cout << "Some required arguments missing" << std::endl;
    std::cout << " Use --help for description of inputs and options"
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
    std::cout << "Error -- number of fluxes requested is zero."
	      << std::endl;
    return 1;
  }
  if (nflux & (nflux-1)) {
    std::cout << "nflux must be power of 2" << std::endl;
    std::cout << " yours is: " << nflux << std::endl;
    return 1;
  }
  if (sigma < 0.0) {
    std::cout << "Invalid noise level: must be >= 0.0" << std::endl;
    return 1;
  }
  if (n0 < 0.0) {
    std::cout << "Invalid (negative) number of sources per area"
	      << std::endl;
    return 1;
  }
  if (nbins == 0) {
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
  if (pixsize >= fwhm/2.0) {
    std::cout << "Insufficient (FWHM/2) beam sampling based on pixel size"
	      << std::endl;
    return 1;
  }
  if (oversample % 2 == 0) {
    std::cout << "Invalid (non-odd) oversampling " << oversample << std::endl;
    return 1;
  }
  if (filterscale < 0.0) {
    std::cout << "Invalid (negative) filter scale " << filterscale << std::endl;
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
	std::cout << "Error reading wisdom file: " << wisdom_file << std::endl;
	return 4;
      }
    }

    if (n0 == 0)
      n0 = model.getBaseN0();
    
    if (verbose) {
      printf("   Beam fwhm:          %0.2f\n", bm.getFWHM());
      printf("   Beam area:          %0.3e\n", bm.getEffectiveArea());
      printf("   Pixel size:         %0.2f\n", pixsize);
      printf("   Flux per area:      %0.2f\n",
	     model.getBaseFluxPerArea());
      printf("   Base N0:            %0.4e\n", model.getBaseN0());
      printf("   N0:                 %0.4e\n", n0);
      printf("   sigma:              %0.4f\n", sigma);
      if (filterscale > 0.0)
	printf("   filter scale:       %0.4f\n", filterscale);
      if (return_log) 
	printf("  Returning log(P(D)) rather than P(D)\n");
      if (oversample != 1)
	printf("oversamp:              %u\n", oversample);
    }

    // Get histogrammed inverse beam
    beamHist inv_bmhist(nbins, filterscale);
    inv_bmhist.fill(bm, nfwhm, pixsize, true, oversample);

    //Get P(D)
    if (verbose) std::cout << "Getting P(D) with transform length: " 
			   << nflux << " and max flux: " 
			   << maxflux << std::endl;
    double base_n0 = model.getBaseN0();
    pfactory.initPD(nflux, sigma, maxflux, base_n0 > n0 ? base_n0 : n0, 
		    model, inv_bmhist);

    if (write_r) {
      if (verbose) std::cout << "Writing R to " << r_file << std::endl;
      pfactory.writeRToFile(r_file);
    }
    
    pfactory.getPD(n0, pd, return_log, true);
    
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
	std::cout << "Error opening output file: " << outputfile
		  << std::endl;
	return 64;
      }
      ofs << pd;
      ofs.close();
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

int getPDDouble(int argc, char **argv) {

  std::string modelfile; //Knot parameters
  double n0; //Number of sources
  double maxflux1, maxflux2; //Maximum fluxes requested
  double fwhm1, fwhm2, nfwhm, pixsize, filterscale; //Calculation params (req)
  double sigma1, sigma2; //Instrument noise
  std::string outputfile; //Ouput pofd option
  unsigned int nflux, nbins, oversample;
  bool has_wisdom, verbose, return_log, write_fits, has_user_pixsize, write_r;
  std::string wisdom_file, r_file;

  //Defaults
  sigma1              = 2e-3;
  sigma2              = 2e-3;
  has_wisdom          = false;
  nflux               = 2048;
  nbins               = 150;
  nfwhm               = 4.5;
  filterscale         = 0.0;
  verbose             = false;
  return_log          = false;
  write_fits          = false;
  has_user_pixsize    = false;
  pixsize             = 3.0;
  write_r             = false;
  oversample          = 1;

  int c;
  int option_index = 0;
  optind = 1; //Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			  &option_index)) != -1) 
    switch(c) {
    case 'f' :
      write_fits = true;
      break;
    case 'F' :
      filterscale = atof(optarg);
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
    case 'o':
      oversample = static_cast<unsigned int>(atoi(optarg));
      break;
    case 'p':
      has_user_pixsize = true;
      pixsize = atof(optarg);
      break;
    case 'r':
      write_r = true;
      r_file = std::string(optarg);
    case '3' :
      sigma1 = atof(optarg);
      break;
    case '4' :
      sigma2 = atof(optarg);
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

  if (optind >= argc - 6) {
    std::cout << "Some required arguments missing" << std::endl;
    std::cout << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  modelfile  = std::string(argv[optind]);
  n0         = atof(argv[optind + 1]);
  fwhm1      = atof(argv[optind + 2]);
  fwhm2      = atof(argv[optind + 3]);
  maxflux1   = atof(argv[optind + 4]);
  maxflux2   = atof(argv[optind + 5]);
  outputfile = std::string(argv[optind + 6]);

  if (!has_user_pixsize)
    pixsize = fwhm1 < fwhm2 ? fwhm1 / 3.0 : fwhm2 / 3.0;

  //Input tests
  if (nflux == 0) {
    std::cout << "Error -- number of fluxes requested is zero."
	      << std::endl;
    return 1;
  }
  if (nflux & (nflux-1)) {
    std::cout << "nflux must be power of 2" << std::endl;
    std::cout << " yours is: " << nflux << std::endl;
    return 1;
  }
  if (sigma1 < 0.0) {
    std::cout << "Invalid noise level (band1): must be >= 0.0" << std::endl;
    return 1;
  }
  if (sigma2 < 0.0) {
    std::cout << "Invalid noise level (band2): must be >= 0.0" << std::endl;
    return 1;
  }
  if (n0 < 0.0) {
    std::cout << "Invalid (negative) number of sources per area"
	      << std::endl;
    return 1;
  }
  if (nbins == 0) {
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
  if (filterscale < 0.0) {
    std::cout << "Invalid (negative) filter scale: " << filterscale
	      << std::endl;
    return 1;
  }
  if (oversample % 2 == 0) {
    std::cout << "Invalid (non-odd) oversampling " << oversample << std::endl;
    return 1;
  }

  try {
    numberCountsDouble model(modelfile);
    doublebeam bm(fwhm1, fwhm2);
    PDFactoryDouble pfactory;
    PDDouble pd;

    if (verbose) pfactory.setVerbose();

    bool success;
    if (has_wisdom) {
      std::cout << "Reading in wisdom file: " << wisdom_file 
		<< std::endl;
      success = pfactory.addWisdom(wisdom_file);
      if (!success) {
	std::cout << "Error reading wisdom file: " << wisdom_file << std::endl;
	return 4;
      }
    }
    
    if (n0 == 0)
      n0 = model.getBaseN0();

    if (verbose) {
      std::pair<double, double> dpr;
      dpr = bm.getFWHM();
      printf("   Beam fwhm1:         %0.2f\n", dpr.first);
      printf("   Beam fwhm2:         %0.2f\n", dpr.second);
      dpr = bm.getEffectiveArea();
      printf("   Beam area1:         %0.3e\n", dpr.first);
      printf("   Beam area2:         %0.3e\n", dpr.second);
      printf("   Pixel size:         %0.2f\n", pixsize);
      printf("   Flux per area1:     %0.2f\n",
	     model.getBaseFluxPerArea1());
      printf("   Flux per area2:     %0.2f\n",
	     model.getBaseFluxPerArea2());
      printf("   Base N0:            %0.4e\n", model.getBaseN0());
      printf("   N0:                 %0.4e\n", n0);
      printf("   sigma1:             %0.4f\n", sigma1);
      printf("   sigma2:             %0.4f\n", sigma2);
      if (filterscale > 0.0)
	printf("   filter scale:       %0.4f\n", filterscale);
      if (return_log) 
	printf("  Returning log(P(D)) rather than P(D)\n");
      if (oversample != 1)
	printf("oversamp:              %u\n", oversample);
    }

    // Get histogrammed inverse beam
    doublebeamHist inv_bmhist(nbins, filterscale);
    inv_bmhist.fill(bm, nfwhm, pixsize, true, oversample);

    //Get P(D)
    if (verbose) std::cout << "Getting P(D) with transform length: " 
			   << nflux << " and max fluxes: " 
			   << maxflux1 << " " << maxflux2 << std::endl;
    double base_n0 = model.getBaseN0();
    pfactory.initPD(nflux, sigma1, sigma2, maxflux1, maxflux2, 
		    base_n0 > n0 ? base_n0 : n0, model, inv_bmhist);

    if (write_r) {
      if (verbose) std::cout << "Writing R to " << r_file << std::endl;
      pfactory.writeRToFile(r_file);
    }

    pfactory.getPD(n0, pd, return_log, true);
    
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
	std::cout << "Error opening output file: " << outputfile
		  << std::endl;
	return 64;
      }
      ofs << pd;
      ofs.close();
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
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'h' :
      std::cout << "NAME" << std::endl;
      std::cout << "\tpofd_coverage_getPD -- get the P(D) for a spline"
		<< std::endl;
      std::cout << "\t type model with a Gaussian beam (1D) or the same type"
		<< " of model" << std::endl;
      std::cout << "\t paired with a log-normal color model in flux_2/flux_1."
		<< std::endl;
      std::cout << std::endl;
      std::cout << "SYNOPSIS" << std::endl;
      std::cout << "\t One-dimensional case:" << std::endl;
      std::cout << "\t  pofd_coverage_getPD [options] modelfile n0 fwhm maxflux"
		<< " outfile" << std::endl; 
      std::cout << std::endl;
      std::cout << "\t Two-dimensional case:" << std::endl;
      std::cout << "\t  pofd_coverage_getPD [options] -d modelfile n0 fwhm1"
		<< " fwhm2" << std::endl;
      std::cout << "\t    maxflux1 maxflux2 outfile" << std::endl; 
      std::cout << std::endl;
      std::cout << "DESCRIPTION" << std::endl;
      std::cout << "\tEvaluates P(D) for the specified model and writes it to" 
		<< std::endl;
      std::cout << "\toutfile.  The number counts model in 1D is a spline"
		<< std::endl;
      std::cout << "\tmodel specified by modelfile, and by the number of"
		<< std::endl;
      std::cout << "\tsources per unit area n0." << std::endl;
      std::cout << std::endl;
      std::cout << "\tIn the 2D case the model is the 1D model in band 1 times"
		<< " a" << std::endl;
      std::cout << "\tLog-Normal distribution in flux2/flux1.  The mu and sigma"
		<< " Log-Normal" << std::endl;
      std::cout << "\tmodel parameters are stored as splines as a function of"
		<< " the" << std::endl;
      std::cout << "\tflux in the first band." << std::endl;
      std::cout << std::endl;
      std::cout << "\tmodelfile should be a text file.  For 1D it consists of"
		<< " nknots" << std::endl;
      std::cout << "\tlines of the form:" << std::endl << std::endl;
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
      std::cout << "\tfwhm is the beam FWHM in arcsec.  The beam is assumed "
		<< "Gaussian. " << std::endl;
      std::cout << "\tIn the 2D case, fwhm1 and fwhm2 are the values for each"
		<< " band." << std::endl;
      std::cout << std::endl;
      std::cout << "\tmaxflux is the maximum flux density generated.  In"
		<< " general" << std::endl;
      std::cout << "\tthe maxflux values will not quite be realized.  Again,"
		<< std::endl;
      std::cout << "\tin the 2D case maxflux1 and maxflux2 are the values in"
		<< " each" << std::endl;
      std::cout << "\tband." << std::endl;
      std::cout << std::endl;
      std::cout << "OPTIONS" << std::endl;
      std::cout << "\t-d, --twod" << std::endl;
      std::cout << "\t\tIf set, the two-dimensional model is used."
		<< std::endl;
      std::cout << "\t-h --help" << std::endl;
      std::cout << "\t\tPrint this message and exit." << std::endl;
      std::cout << "\t-f, --fits" << std::endl;
      std::cout << "\t\tWrite output as a fits file rather than text."
		<< std::endl;
      std::cout << "\t-F, --filtscale VALUE" << std::endl;
      std::cout << "\t\tRadius of high-pass filter in arcseconds. If zero,"
		<< std::endl;
      std::cout << "\t\tno filtering is applied (def: 0)." << std::endl;
      std::cout << "\t-l, --log" << std::endl;
      std::cout << "\t\tReturn the log P(D) rather than the P(D)."
		<< std::endl;
      std::cout << "\t-n, --nflux value" << std::endl;
      std::cout << "\t\tThe number of requested fluxes along each dimension."
		<< std::endl;
      std::cout << "\t\tAlso sets the transform size used. (def: 131072 in 1D,)"
		<< std::endl;
      std::cout << "\t\tand 2048 in 2D)." << std::endl;
      std::cout << "\t-N, --nfwhm value" << std::endl;
      std::cout << "\t\tNumber of beam FWHM out to go when computing beam."
		<< "(def: 3.5)" << std::endl;
      std::cout << "\t--nbins value" << std::endl;
      std::cout << "\t\tNumber of bins to use in histogrammed beam. (def: 80)"
		<< std::endl;
      std::cout << "\t-o, --oversample VALUE" << std::endl;
      std::cout << "\t\tAmount to oversample the beam; must be odd integer."
		<< " (def: 1)" << std::endl;
      std::cout << "\t-p, --pixsize value" << std::endl;
      std::cout << "\t\tPixel size in arcsec. (def: FWHM/3.0)" << std::endl;
      std::cout << "\t-r, --rfile FILENAME" << std::endl;
      std::cout << "\t\tWrite the R used to this file as text." << std::endl;
      std::cout << "\t-v, --verbose" << std::endl;
      std::cout << "\t\tPrint informational messages while running"
		<< std::endl;
      std::cout << "\t-V, --version" << std::endl;
      std::cout << "\t\tOutput version number and exit" << std::endl;
      std::cout << "\t-w, --wisdom wisdomfile" << std::endl;
      std::cout << "\t\tName of wisdom file (prepared with fftw-wisdom)." 
		<< std::endl;
      std::cout << "\tONE-D MODEL OPTIONS" << std::endl;
      std::cout << "\t-s, --sigma VALUE" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise (def: 0.002)" << std::endl;
      std::cout << "\tTWO-D MODEL OPTIONS" << std::endl;
      std::cout << "\t--sigma1 NOISE" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise, band 1 (def: 0.002)." 
		<< std::endl;
      std::cout << "\t--sigma2 NOISE" << std::endl;
      std::cout << "\t\tThe assumed per-pixel noise, band 2 (def: 0.002)." 
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
    return getPDSingle(argc, argv);
  else
    return getPDDouble(argc, argv);
}
