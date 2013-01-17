#include<iostream>
#include<string>
#include<iomanip>
#include<vector>
#include<fstream>
#include<cstdlib>

#include<getopt.h>

#include<simManager.h>
#include<pofdExcept.h>
#include<global_settings.h>

static struct option long_options[] = {
  {"help", no_argument, 0, 'h'},
  {"bin",no_argument,0,'b'},
  {"double", no_argument, 0, 'd'},
  {"esmooth",required_argument,0,'e'},
  {"fftsize",required_argument,0,'f'},
  {"nbins",required_argument,0,'1'},
  {"nsims",required_argument,0,'n'},
  {"nolike",no_argument,0,'N'},
  {"nlike",required_argument,0,'2'},
  {"n0initrange",required_argument,0,'3'},
  {"n0rangefrac",required_argument,0,'4'},
  {"oversample",required_argument,0,'o'},
  {"sigma",required_argument,0,'s'},
  {"seed",required_argument,0,'S'},
  {"verbose",no_argument,0,'v'},
  {"version",no_argument,0,'V'},
  {"wisdom",required_argument,0,'w'},
  {0,0,0,0}
};

char optstring[] = "hbde:f:1:n:N2:3:4:o:s:S:vVw:";

int runSimSingle(int argc, char **argv) {

  std::string modelfile; //Base model file
  double n0; //Model params
  double fwhm; //Calculation params (req)
  double esmooth; //Extra smoothing amount
  unsigned int nsims, nlike, n1, n2, fftsize, nbins, oversample;
  double pixsize, n0rangefrac, n0initrange;
  double sigma; //Instrument noise, unsmoothed
  std::string outputfile; //Ouput pofd option
  bool verbose, has_wisdom, has_user_seed, use_binning, map_like;
  unsigned long long int seed;
  std::string wisdom_file;

  //Defaults
  fftsize             = 131072;
  nsims               = 100;
  nlike               = 500;
  n0rangefrac         = 0.1;
  n0initrange         = 0.15;
  sigma               = 0.002;
  verbose             = false;
  has_wisdom          = false;
  esmooth             = 0.0;
  has_user_seed       = false;
  seed                = 1024;
  oversample          = 1;
  nbins               = 1000;
  use_binning         = false;
  map_like            = true;
  
  int c;
  int option_index = 0;
  optind = 1; //Rewind 
  while ( ( c = getopt_long(argc,argv,optstring,long_options,
			    &option_index ) ) != -1 ) 
    switch(c) {
    case 'b' :
      use_binning = true;
      break;
    case 'e' :
      esmooth = atof(optarg);
      break;
    case 'f' :
      fftsize = atoi(optarg);
      break;
    case '1' :
      nbins = atoi(optarg);
      break;
    case 'n' :
      nsims = atoi(optarg);
      break;
    case 'N' :
      map_like = false;
      break;
    case '2' :
      nlike = atoi(optarg);
      break;
    case '3' :
      n0initrange = atof(optarg);
      break;
    case '4' :
      n0rangefrac = atof(optarg);
      break;
    case 'o' :
      oversample = atoi(optarg);
      break;
    case 's' :
      sigma = atof(optarg);
      break;
    case 'S' :
      has_user_seed = true;
      seed = static_cast<unsigned long long int>( atoi(optarg) );
      break;
    case 'v' :
      verbose = true;
      break;
    case 'w' :
      has_wisdom = true;
      wisdom_file = std::string( optarg );
      break;
    }

  if (optind >= argc-6 ) {
    std::cerr << "Some required arguments missing" << std::endl;
    std::cerr << " Use --help for description of inputs and options"
	      << std::endl;
    return 1;
  }
  n0         = atof(argv[optind]);
  modelfile  = std::string(argv[optind + 1]);
  fwhm       = atof(argv[optind + 2]);
  pixsize    = atof(argv[optind + 3]);
  n1         = atoi(argv[optind + 4]);
  n2         = atoi(argv[optind + 5]);
  outputfile = std::string(argv[optind + 6]);

  if (sigma < 0.0) {
    std::cerr << "Invalid instrument noise level: must be >= 0.0" << std::endl;
    return 1;
  }
  if (n0initrange <= 0.0) {
    std::cerr << "Invalid n0initrange: must be > 0" << std::endl;
    return 1;
  }
  if (n0initrange >= 1.0) {
    std::cerr << "Invalid n0initrange: must be < 1" << std::endl;
    return 1;
  }
  if (fwhm <= 0.0) {
    std::cerr << "Invalid (non-positive) FWHM" << std::endl;
    return 1;
  }
  if (esmooth < 0.0) {
    std::cerr << "Invalid (negative) additional smoothing" << std::endl;
    return 1;
  }
  if (use_binning && nbins == 0) {
    std::cerr << "Invalid (zero) number of bins" << std::endl;
    return 1;
  }
  if (oversample == 0) {
    std::cerr << "Invalid (non-positive) oversampling" << std::endl;
    return 1;
  }
  if (n0 <= 0.0) {
    std::cerr << "Invalid (non-positve) n0" << std::endl;
    return 1;
  }
  if (pixsize <= 0.0) {
    std::cerr << "Invalid (non-positve) pixsize" << std::endl;
    return 1;
  }
  if (n1*n2 == 0) {
    std::cerr << "Simulated image has zero size" << std::endl;
    return 1;
  }
  if (map_like && n0rangefrac <= 0.0) {
    std::cerr << "Invalid n0rangefrac: must be > 0" << std::endl;
    return 1;
  }
  if (map_like && n0rangefrac >= 1.0) {
    std::cerr << "Invalid n0rangefrac: must be < 1" << std::endl;
    return 1;
  }

  try {
    if (verbose) {
      double area = n1*n2*std::pow(pixsize/3600.0,2);
      printf("   base model file:    %s\n", modelfile.c_str());
      printf("   nsims:              %u\n",nsims);
      printf("   n0initrange         %0.3f\n",n0initrange);
      printf("   Beam fwhm:          %0.2f\n",fwhm);
      printf("   Pixel size:         %0.1f\n",pixsize);
      printf("   Number of pixels:   %u x %u\n",n1,n2);
      printf("   Area:               %0.2f\n",area);
      printf("   N0:                 %0.3e\n",n0);
      printf("   sigma:              %0.4f\n",sigma);
      printf("   fftsize:            %u\n",fftsize);
      if (esmooth > 0) 
	printf("   esmooth:            %0.2f\n",esmooth);
      if (oversample > 1)
	printf("   oversampling:       %u\n",oversample);
      if (use_binning)
	printf("   nbins:              %u\n",nbins);
      if (map_like) {
	printf("   nlike:              %u\n",nlike);
	printf("   n0rangefrac         %0.3f\n",n0rangefrac);
      }
    }

    simManager sim(modelfile, nsims, n0initrange, map_like, nlike, 
		   n0rangefrac, fftsize, n1, n2, pixsize, fwhm, sigma, n0, 
		   esmooth, oversample, use_binning, nbins);
    if (has_wisdom) sim.addWisdom(wisdom_file);
    if (has_user_seed) sim.setSeed(seed);

    sim.doSims(verbose);

    //Write it
    if (verbose) std::cout << "Writing simulation results to " << outputfile 
			   << std::endl;
    int status = sim.writeToFits(outputfile);
    if (status != 0) return status;
  } catch (const pofdExcept& ex) {
    std::cerr << "Error encountered" << std::endl;
    std::cerr << ex << std::endl;
    return 8;
  } catch (const std::bad_alloc& ba) {
    std::cerr << "Bad allocation error: " << ba.what() << std::endl;
    return 16;
  } 

  return 0;
}


////////////////////////////////////////


int main(int argc, char **argv) {
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
      std::cerr << "\tpofd_coverage_runSim -- make a set of "
		<< "simulated images for" << std::endl;
      std::cerr << "\t a broken power law model with a"
		<< " Gaussian beam and" << std::endl;
      std::cerr << "\t measure the number of objects per sq deg from them."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "SYNOPSIS" << std::endl;
      std::cerr << "\t pofd_coverage_runSim [options] n0 modelfile fwhm pixsize "
		<< "n1 n2" << std::endl;
      std::cerr << "\t  outputfile" << std::endl;
      std::cerr << std::endl;
      std::cerr << "DESCRIPTION" << std::endl;
      std::cerr << "\tCreates simulated images for a given model, computes the"
		<< " P(D)" << std::endl;
      std::cerr << "\tfor a range of input n0 values, and finds the best"
		<< " fitting value" << std::endl;
      std::cerr << "\tof n0 for each image.  Optionally, this is then followed"
		<< " by" << std::endl;
      std::cerr << "\tmapping out a range of n0 values around the best fit and"
		<< " storing" << std::endl;
      std::cerr << "\tthe resulting log likelihood.  The results are written to"
		<< " outputfile." << std::endl;
      std::cerr << std::endl;
      std::cerr << "\tThe number counts model is a broken power law model"
		<< std::endl;
      std::cerr << "\tspecified by modelfile and by the number of sources"
		<< std::endl;
      std::cerr << "\tper unit area n0." << std::endl;
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
      std::cerr << "\tfwhm is the beam FWHM in arcsec , and the beams is "
		<< "Gaussian."
		<< std::endl;
      std::cerr << std::endl;
      std::cerr << "\tPixsize gives the pixel size (in arcsec), n1 and n2 give"
		<< " the" << std::endl;
      std::cerr << "\tnumber of pixels along each dimension in the simulated "
		<< "data." << std::endl;
      std::cerr << std::endl;
      std::cerr << "OPTIONS" << std::endl;
      std::cerr << "\t-h --help" << std::endl;
      std::cerr << "\t\tPrint this message and exit." << std::endl;
      std::cerr << "\t-b, --bin" << std::endl;
      std::cerr << "\t\tBin the simulated image for the likelihood calculation."
		<< std::endl;
      std::cerr << "\t-f, --fftsize FFTSIZE" << std::endl;
      std::cerr << "\t\tSize of FFT to use when computing P(D) (def: 131072 in"
		<< std::endl;
      std::cerr << "\t\tone dimension and 2048 in two dimensions." 
		<< std::endl;
      std::cerr << "\t--nbins NBINS" << std::endl;
      std::cerr << "\t\tNumber of bins to use if binning simulated image."
		<< std::endl;
      std::cerr << "\t-n, --nsims NSIMS" << std::endl;
      std::cerr << "\t\tThe number of simulations to do (def: 100)." 
		<< std::endl;
      std::cerr << "\t-N, --nolike" << std::endl;
      std::cerr << "\t\tDo not map out the likelihood values." << std::endl;
      std::cerr << "\t--nlike NLIKE" << std::endl;
      std::cerr << "\t\tThe number of likelihoods to compute for each sim if" 
		<< std::endl;
      std::cerr << "\t\t --nolike is not set (def: 500)." << std::endl;
      std::cerr << "\t--n0initrange N0INITRANGE" << std::endl;
      std::cerr << "\t\tFractional range used to bracket likelihood in "
		<< "initial" << std::endl;
      std::cerr << "\t\tsearch for maximum likelihood (def: 0.15)."
		<< std::endl;
      std::cerr << "\t--n0rangefrac N0RANGEFRAC" << std::endl;
      std::cerr << "\t\tThe fractional range in n0 to explore in each "
		<< "direction" << std::endl;
      std::cerr << "\t\tif doing likelihood map (def: 0.1)."
		<< std::endl;
      std::cerr << "\t-S, --seed SEED" << std::endl;
      std::cerr << "\t\tSet user specified seed, otherwise taken from time."
		<< std::endl;
      std::cerr << "\t-v, --verbose" << std::endl;
      std::cerr << "\t\tPrint informational messages while running"
		<< std::endl;
      std::cerr << "\t-V, --version" << std::endl;
      std::cerr << "\t\tOutput version number and exit" << std::endl;
      std::cerr << "\t-w, --wisdom wisdomfile" << std::endl;
      std::cerr << "\t\tName of wisdom file (prepared with fftw-wisdom)." 
		<< std::endl;
      std::cerr << "ONE-DIMENSIONAL OPTIONS" << std::endl;
      std::cerr << "\t-e, --esmooth ESMOOTH" << std::endl;
      std::cerr << "\t\tExtra smoothing FWHM (in arcsec)" << std::endl;
      std::cerr << "\t-o, --oversample VALUE" << std::endl;
      std::cerr << "\t\tAmount of oversampling to use (integral) when " 
		<< "generating" << std::endl;
      std::cerr << "\t\timage.  The data is then down-binned to the specified"
		<< "size." << std::endl;
      std::cerr << "\t\tThe default is to apply no oversampling." << std::endl;
      std::cerr << "\t-s, --sigma noise" << std::endl;
      std::cerr << "\t\tThe assumed per-pixel noise (assumed Gaussian).  This"
		<< " is" << std::endl;
      std::cerr << "\t\tthe level before any extra smoothing (def: 0.002)."
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
    return runSimSingle(argc,argv);
  else {
    std::cerr << "2D model not supported" << std::endl;
    return 1;
  }
}