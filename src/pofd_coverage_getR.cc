#include<iostream>

#include<getopt.h>

#include "../include/global_settings.h"
#include "../include/utility.h"
#include "../include/beam.h"
#include "../include/doublebeam.h"
#include "../include/numberCounts.h"
#include "../include/numberCountsDouble.h"
#include "../include/pofdExcept.h"

//All sub-parses have to have the same options to avoid
// getting warnings.  So we give them all the same long_options,
// but then only process the appropriate ones, ignoring the rest
//Yes, this is a bit complicated and error prone, but such is life
static struct option long_options[] = {
  {"help", no_argument, 0, 'h'},
  {"double", no_argument, 0, 'd'},
  {"filterscale", required_argument, 0, 'F'},
  {"nbins", required_argument, 0, 'n'},
  {"nfwhm", required_argument, 0, 'N'},
  {"oversamp", required_argument, 0, 'o'},
  {"verbose", no_argument, 0, 'v'},
  {"version", no_argument, 0, 'V'}, 
  {0,0,0,0}
};

char optstring[] = "hdF:n:N:o:vV";

//One-D version
int getRSingle(int argc, char** argv) {

  std::string modelfile; //Init file (having model we want)
  std::string outfile; //File to write to
  bool verbose;
  double n0, nfwhm, pixsize, maxflux, filterscale, fwhm;
  unsigned int nflux, nbins, oversamp;

  // Defaults
  verbose = false;
  nfwhm = 40.0;
  nbins = 120;
  filterscale = 0.0;
  oversamp = 1;

  int c;
  int option_index = 0;
  optind = 1; //!< Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			    &option_index)) != -1) 
    switch(c) {
    case 'F':
      filterscale = atof(optarg);
      break;
    case 'n':
      nbins = atoi(optarg);
      break;
    case 'N':
      nfwhm = atof(optarg);
      break;
    case 'o':
      oversamp = atoi(optarg);
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc - 6) {
    std::cerr << "Required arguments missing" << std::endl;
    return 1;
  }
  modelfile = std::string(argv[optind]);
  n0 = atof(argv[optind+1]);
  fwhm = atof(argv[optind+2]);
  pixsize = atof(argv[optind+3]);
  maxflux = atof(argv[optind+4]);
  nflux = static_cast<unsigned int>(atoi(argv[optind+5]));
  outfile = std::string(argv[optind+6]);

  if (nflux == 0) {
    std::cout << "Error -- number of fluxes requested is zero."
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
  if (filterscale < 0.0) {
    std::cout << "Invalid (negative) filter scale " << filterscale << std::endl;
    return 1;
  }
  if (oversamp % 2 == 0) {
    std::cout << "Invalid (non-odd) oversampling " << oversamp << std::endl;
    return 1;
  }
  if (maxflux <= 0.0) {
    std::cout << "Invalid (non-positive) maxflux " << maxflux << std::endl;
    return 1;
  }
  if (n0 < 0.0) {
    std::cout << "Invalid (non-positive) n0 " << n0 << std::endl;
    return 1;
  }

  double *R = NULL;
  try {
    numberCounts model(modelfile);
    beam bm(fwhm);
    beamHist inv_bmhist(nbins, filterscale);
    inv_bmhist.fill(bm, nfwhm, pixsize, true, oversamp);
    
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
      if (filterscale > 0.0)
	printf("   filter scale:       %0.4f\n", filterscale);
      if (oversamp != 1)
	printf("   oversamp:           %u\n", oversamp);
    }

    // Get R
    R = new double[nflux];
    model.getR(nflux, 0, maxflux, inv_bmhist, R);
    
    // Adjust for N0
    double n0fac = n0 / model.getBaseN0();
    if (n0fac != 1)
      for (unsigned int i = 0; i < nflux; ++i)
	R[i] *= n0fac;

    // Write
    double dflux;
    if (nflux > 1)
      dflux = maxflux / static_cast<double>(nflux - 1);
    else 
      dflux = 0.0;
    FILE *fp;
    fp = fopen( outfile.c_str(),"w");
    if (!fp) {
      std::cerr << "Failed to open output file" << std::endl;
      return 128;
    }
    fprintf(fp, "#%-11s   %-12s\n", "Flux", "R");
    for (unsigned int i = 0; i < nflux; ++i) 
      fprintf(fp, "%12.6e   %15.9e\n", i*dflux, R[i]);
    fclose(fp);
    
    delete[] R; 
  } catch (const pofdExcept& ex) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (R != NULL) delete[] R;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    if (R != NULL) delete[] R;
    return 16;
  }
  return 0;
}

//One-D version
int getRDouble(int argc, char** argv) {

  std::string modelfile; //Init file (having model we want)
  std::string outfile; //File to write to
  bool verbose;
  double n0, nfwhm, pixsize, maxflux1, maxflux2, filterscale, fwhm1, fwhm2;
  unsigned int nflux, nbins, oversamp;

  // Defaults
  verbose = false;
  nfwhm = 40.0;
  nbins = 150;
  filterscale = 0.0;
  oversamp = 1;

  int c;
  int option_index = 0;
  optind = 1; //!< Reset parse
  while ((c = getopt_long(argc,argv,optstring,long_options,
			    &option_index)) != -1) 
    switch(c) {
    case 'F':
      filterscale = atof(optarg);
      break;
    case 'n':
      nbins = atoi(optarg);
      break;
    case 'N':
      nfwhm = atof(optarg);
      break;
    case 'o':
      oversamp = atoi(optarg);
      break;
    case 'v' :
      verbose = true;
      break;
    }

  if (optind >= argc - 8) {
    std::cerr << "Required arguments missing" << std::endl;
    return 1;
  }
  modelfile = std::string(argv[optind]);
  n0 = atof(argv[optind+1]);
  fwhm1 = atof(argv[optind+2]);
  fwhm2 = atof(argv[optind+3]);
  pixsize = atof(argv[optind+4]);
  maxflux1 = atof(argv[optind+5]);
  maxflux2 = atof(argv[optind+6]);
  nflux = static_cast<unsigned int>(atoi(argv[optind+7]));
  outfile = std::string(argv[optind+8]);

  if (nflux == 0) {
    std::cout << "Error -- number of fluxes requested is zero."
	      << std::endl;
    return 1;
  }
  if (nbins == 0) {
    std::cout << "Invalid (non-positive) number of beam histogram bins"
	      << std::endl;
    return 1;
  }
  if (fwhm1 <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM1" << fwhm1 << std::endl;
    return 1;
  }
  if (fwhm2 <= 0.0) {
    std::cout << "Invalid (non-positive) FWHM2" << fwhm2 << std::endl;
    return 1;
  }
  if (nfwhm <= 0) {
    std::cout << "Invalid (non-positive) number of beam FWHMs"
	      << std::endl;
    return 1;
  }
  if (filterscale < 0.0) {
    std::cout << "Invalid (negative) filter scale " << filterscale << std::endl;
    return 1;
  }
  if (oversamp % 2 == 0) {
    std::cout << "Invalid (non-odd) oversampling " << oversamp << std::endl;
    return 1;
  }
  if (maxflux1 <= 0.0) {
    std::cout << "Invalid (non-positive) maxflux1 " << maxflux1 << std::endl;
    return 1;
  }
  if (maxflux2 <= 0.0) {
    std::cout << "Invalid (non-positive) maxflux1 " << maxflux2 << std::endl;
    return 1;
  }
  if (n0 < 0.0) {
    std::cout << "Invalid (non-positive) n0 " << n0 << std::endl;
    return 1;
  }

  // Minflux is 0

  double *R = NULL;
  double *flux1 = NULL;
  double *flux2 = NULL;
  try {
    numberCountsDouble model(modelfile);
    doublebeam bm(fwhm1, fwhm2);
    doublebeamHist inv_bmhist(nbins, filterscale);
    inv_bmhist.fill(bm, nfwhm, pixsize, true, oversamp);

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
      if (filterscale > 0.0)
	printf("   filter scale:       %0.4f\n", filterscale);
      if (oversamp != 1)
	printf("   oversamp:           %u\n", oversamp);
    }

    // Set up fluxes
    double dflux1, dflux2;
    if (nflux > 1) {
      dflux1 = maxflux1 / static_cast<double>(nflux - 1);
      dflux2 = maxflux2 / static_cast<double>(nflux - 1);
    } else {
      dflux1 = 0.0;
      dflux2 = 0.0;
    }
    flux1 = new double[nflux];
    for (unsigned int i = 0; i < nflux; ++i)
      flux1[i] = dflux1 * static_cast<double>(i);
    flux2 = new double[nflux];
    for (unsigned int i = 0; i < nflux; ++i)
      flux2[i] = dflux2 * static_cast<double>(i);

    // Get R
    R = new double[nflux * nflux];
    model.getR(nflux, flux1, nflux, flux2, inv_bmhist, R);
    delete[] flux1;
    delete[] flux2;

    // Adjust for N0
    double n0fac = n0 / model.getBaseN0();
    if (n0fac != 1)
      for (unsigned int i = 0; i < nflux * nflux; ++i)
	R[i] *= n0fac;

    // Write

    FILE *fp;
    fp = fopen(outfile.c_str(), "w");
    if (!fp) {
      std::cerr << "Failed to open output file" << std::endl;
      return 128;
    }
    fprintf(fp,"#%4u %4u\n", nflux, nflux);
    fprintf(fp,"#minflux1: %12.6e dflux1: %12.6e\n", 0., dflux1);
    fprintf(fp,"#minflux2: %12.6e dflux2: %12.6e\n", 0., dflux2);
    for (unsigned int i = 0; i < nflux; ++i) {
      for (unsigned int j = 0; j < nflux - 1; ++j)
        fprintf(fp,"%13.7e ",R[nflux * i + j]);
      fprintf(fp,"%13.7e\n",R[nflux * i + nflux - 1]);
    }
    fclose(fp);

    delete[] R; 
  } catch (const pofdExcept& ex) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    if (R != NULL) delete[] R;
    if (flux1 != NULL) delete[] flux1;
    if (flux2 != NULL) delete[] flux2;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    if (R != NULL) delete[] R;
    if (flux1 != NULL) delete[] flux1;
    if (flux2 != NULL) delete[] flux2;
    return 16;
  }
  return 0;
}


////////////////////////////////////////

int main(int argc, char** argv) {
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
      std::cerr << "NAME" << std::endl;
      std::cerr << "\tpofd_coverage_getR -- get R for a number counts model."
		<< "  Both" << std::endl;
      std::cerr << "\tone-dimensional and two-dimensional models are supported."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\tEither" << std::endl;
      std::cerr << std::endl;
      std::cerr << "\t pofd_coverage_getR [options] modelfile n0 fwhm pixsize"
		<< " maxflux" << std::endl;
      std::cerr << "\t\tnflux outfile" << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tfor the 1D case or" << std::endl;
      std::cerr << std::endl;
      std::cerr << "\t pofd_coverage_getR -d [options] modelfile n0 fwhm1 fwhm2 "
		<< "pixsize" << std::endl;
      std::cerr << "\t\tmaxflux1 maxflux2 nflux outfile" << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tfor the 2D case." << std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tEvaluates R for the model in initfile using the P(D)"
		<< " formalism and" << std::endl;
      std::cerr << "\twrites it to outfile.  The 1D model is a log-space "
		<< "spline" << std::endl;
      std::cerr << "\tmodel for the number counts, and the 2D model is the 1D" 
		<< " spline" << std::endl;
      std::cerr << "\tmodel times a log-normal color function for the second"
		<< " band," << std::endl;
      std::cerr << "\twith the log-space variance and mean color stored as"
		<< " splines" << std::endl;
      std::cerr << "\tin the flux of the first band." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tmodelfile is a text file specifying the model; the exact"
		<< " details" << std::endl;
      std::cerr << "\t(given below) depend on whether the 1D or 2D case is"
		<< " being used." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tFor the 1D case, modelfile is a text file giving the "
		<< "positions" << std::endl;
      std::cerr << "\tof the spline knots and their values in the format"
		<< " knotflux value." << std::endl;
      std::cerr << "\tAdditional elements on each line are ignored."
		<< std::endl;
      std::cerr << "\tFor the 2D case, modelfile is a text file giving the "
		<< "positions" << std::endl;
      std::cerr << "\tof the knot points and their values, followed by the "
		<< "sigma" << std::endl;
      std::cerr << "\tknot positions and their values, then likewise for the "
		<< "colour" << std::endl;
      std::cerr << "\toffset.  The format is three numbers on the first line, "
		<< "giving" << std::endl;
      std::cerr << "\tthe number of number count knots, sigma knots, and "
		<< "offset knots," << std::endl;
      std::cerr << "\tfollowed by a number of lines again with the format"
		<< std::endl;
      std::cerr << "\tknotpos value.  The sigmas and offsets are in log space."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tn0 is the number of sources per square degree.  If set"
		<< " to" << std::endl;
      std::cerr << "\tzero, then the value from the input model is used "
		<< "directly." << std::endl;
      std::cerr << "\tOtherwise, the model is scaled to match this value."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tfwhm is the beam FWHM (in arcsec) in the 1D case, and"
		<< std::endl;
      std::cerr << "\tfwhm1, fwhm2 are the FWHM values for each band in the 2D"
		<< "case." << std::endl;
      std::cerr << "\tpixscale is the pixel scale (in arcsec)." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tminflux, maxflux give the minimum and maximum"
		<< " flux density" << std::endl;
      std::cerr << "\tto evaluate R for in the 1D case.  For the 2D case this" 
		<< std::endl;
      std::cerr << "\tis extended to the minimum and maximum in each band."
		<< " nflux" << std::endl;
      std::cerr << "\tis the number of fluxes to generate; in the 2D case along"
		<< " each" << std::endl;
      std::cerr << "\tdimension." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tIn both cases the output R is written to outfile as"
		<< " text." << std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-d, --double" << std::endl;
      std::cerr << "\t\tUse the 2D model." << std::endl;
      std::cerr << "\t-F, --filterscale VALUE" << std::endl;
      std::cerr << "\t\tHigh-pass filter scale, in arcsec.  Zero means no"
		<< " filtering" << std::endl;
      std::cerr << "\t\tis applied (def: 0)" << std::endl;
      std::cerr << "\t-n, --nbins VALUE" << std::endl;
      std::cerr << "\t\tNumber of beam histogram bins (def: 120)" << std::endl;
      std::cerr << "\t-N, --nfwhm VALUE" << std::endl;
      std::cerr << "\t\tNumber of FWHM to go out in beam representation. "
		<< "(def: 40.0)" << std::endl;
      std::cerr << "\t-o, --oversamp VALUE" << std::endl;
      std::cerr << "\t\tOversampling of pixels used to generate beam. One means"
		<< std::endl;
      std::cerr << "\t\tno oversampling.  Must be odd (def: 1)" << std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput the version number and exit." << std::endl;
      std::cerr << std::endl;
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

  if (! twod)
    return getRSingle(argc,argv);
  else
    return getRDouble(argc,argv);
}
