//Model testing
#include<string>
#include<iostream>
#include<cmath>

#include<gtest/gtest.h>

#include "../include/global_settings.h"
#include "../include/pofdExcept.h"
#include "../include/beam.h"
#include "../include/doublebeam.h"
#include "../include/numberCounts.h"
#include "../include/numberCountsDouble.h"

////////////////////////////////////
// beam
TEST(beam1DTest, Basic) {
  beam bm(15.0);
  
  EXPECT_NEAR(15.0, bm.getFWHM(), 1e-5) << "Unexpected beam FWHM";
  EXPECT_NEAR(159701.110401, bm.getRhoSq(), 1e-3) << "Unexpected rho2";
  EXPECT_NEAR(1.96717e-5, bm.getEffectiveArea(), 1e-8) << 
    "Unexpected beam area";
  EXPECT_NEAR(9.83585100e-6, bm.getEffectiveAreaSq(), 1e-8) << 
    "Unexpected beam area";

  bm.setFWHM(22.1);
  EXPECT_NEAR(22.1, bm.getFWHM(), 1e-5) << "Unexpected beam FWHM after setFWHM";
  EXPECT_NEAR(73570.872505, bm.getRhoSq(), 1e-3) << "Unexpected rho2";
  EXPECT_NEAR(4.2701582e-5, bm.getEffectiveArea(), 1e-7) << 
    "Unexpected beam area";
  EXPECT_NEAR(2.1350791e-5, bm.getEffectiveAreaSq(), 1e-8) << 
    "Unexpected beam area";
}

TEST(beam1DTest, BeamFactor) {
  const double fwhm = 15;
  const double pixsize = 5.0;
  beam bm(fwhm);

  const unsigned int n = 5;
  double bmfac[n];

  EXPECT_THROW(bm.getBeamFac(0, pixsize, bmfac), pofdExcept) <<
    "Asking getBeamFac for no pixels should throw exception";
  EXPECT_THROW(bm.getBeamFac(4, pixsize, bmfac), pofdExcept) <<
    "Asking getBeamFac for even number of pixels should throw exception";
  EXPECT_THROW(bm.getBeamFac(n, -pixsize, bmfac), pofdExcept) <<
    "Asking getBeamFac for negative pixel size should throw exception";
  EXPECT_THROW(bm.getBeamFac(n, pixsize, NULL), pofdExcept) <<
    "Null factor array should throw exception";

  bm.getBeamFac(n, pixsize, bmfac);
  const double expfac[n] = {0.29163225989, 0.734867246, 1.0, 
			    0.734867246, 0.29163225989};
  for (unsigned int i = 0; i < n; ++i)
    EXPECT_NEAR(expfac[i], bmfac[i], 1e-4) <<
      "Got unexpected beam factor at pixel number " << i << " for " <<
      pixsize << " arcsec pix";

  //Try with different pixel size
  bm.getBeamFac(n, 3.0, bmfac);
  const double expfac2[n] = {0.6417129487, 0.8950250709, 1.0, 
			     0.8950250709, 0.6417129487};
  for (unsigned int i = 0; i < n; ++i)
    EXPECT_NEAR(expfac2[i], bmfac[i], 1e-4) <<
      "Got unexpected beam factor at pixel number " << i << " for 3 arcsec pix";

}

TEST(beam1DTest, GetBeam) {
  const double fwhm = 15;
  const double pixsize = 5.0;
  beam bm(fwhm);

  const unsigned int n = 5;
  double bmarr[n * n];

  EXPECT_THROW(bm.getBeam(0, pixsize, bmarr), pofdExcept) <<
    "Asking getBeam for no pixels should throw exception";
  EXPECT_THROW(bm.getBeam(4, pixsize, bmarr), pofdExcept) <<
    "Asking getBeam for even number of pixels should throw exception";
  EXPECT_THROW(bm.getBeam(n, -pixsize, bmarr), pofdExcept) <<
    "Asking getBeam for negative pixel size should throw exception";
  EXPECT_THROW(bm.getBeam(n, pixsize, NULL), pofdExcept) <<
    "Null factor array should throw exception";

  bm.getBeam(n, pixsize, bmarr);
  const double expfac[n] = {0.29163225989, 0.734867246, 1.0, 
			    0.734867246, 0.29163225989};
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
    EXPECT_NEAR(expfac[i] * expfac[j], bmarr[i * n + j], 1e-4) <<
      "Got unexpected beam value at pixel number " << i << " " << j;
}


////////////////////////////////////
// doublebeam
TEST(beam2DTest, Basic) {
  doublebeam bm(15.0, 20.0);
  
  EXPECT_NEAR(15.0, bm.getFWHM1(), 1e-5) << "Unexpected beam FWHM1";
  EXPECT_NEAR(159701.110401, bm.getRhoSq1(), 1e-3) << "Unexpected rho2 band 1";
  EXPECT_NEAR(1.96717e-5, bm.getEffectiveArea1(), 1e-8) << 
    "Unexpected beam area band 1";
  EXPECT_NEAR(9.83585100e-6, bm.getEffectiveAreaSq1(), 1e-8) << 
    "Unexpected beam area band 1";

  EXPECT_NEAR(20.0, bm.getFWHM2(), 1e-5) << "Unexpected beam FWHM2";
  EXPECT_NEAR(89831.874601, bm.getRhoSq2(), 1e-3) << "Unexpected rho2 band 2";
  EXPECT_NEAR(3.49719147e-5, bm.getEffectiveArea2(), 1e-8) << 
    "Unexpected beam area band 2";
  EXPECT_NEAR(1.7485957337e-5, bm.getEffectiveAreaSq2(), 1e-8) << 
    "Unexpected beam area band 2";


  bm.setFWHM(19.0, 22.0);
  EXPECT_NEAR(19.0, bm.getFWHM1(), 1e-5) << 
    "Unexpected beam FWHM band 1 after setFWHM";
  EXPECT_NEAR(99536.7031585, bm.getRhoSq1(), 1e-3) << 
    "Unexpected rho2 band 1 after setFWHM";
  EXPECT_NEAR(22.0, bm.getFWHM2(), 1e-5) << 
    "Unexpected beam FWHM band 2 after setFWHM";
  EXPECT_NEAR(74241.2186782, bm.getRhoSq2(), 1e-3) << 
    "Unexpected rho2 band 2 after setFWHM";
}

TEST(beam2DTest, BeamFactor) {
  const double fwhm1 = 15;
  const double fwhm2 = 20;
  const double pixsize = 5.0;
  doublebeam bm(fwhm1, fwhm2);

  const unsigned int n = 5;
  double bmfac[n];

  EXPECT_THROW(bm.getBeamFac1(0, pixsize, bmfac), pofdExcept) <<
    "Asking getBeamFac1 for no pixels should throw exception";
  EXPECT_THROW(bm.getBeamFac1(4, pixsize, bmfac), pofdExcept) <<
    "Asking getBeamFac1 for even number of pixels should throw exception";
  EXPECT_THROW(bm.getBeamFac1(n, -pixsize, bmfac), pofdExcept) <<
    "Asking getBeamFac1 for negative pixel size should throw exception";
  EXPECT_THROW(bm.getBeamFac1(n, pixsize, NULL), pofdExcept) <<
    "Null factor array should throw exception";
  bm.getBeamFac1(n, pixsize, bmfac);
  const double expfac[n] = {0.29163225989, 0.734867246, 1.0, 
			    0.734867246, 0.29163225989};
  for (unsigned int i = 0; i < n; ++i)
    EXPECT_NEAR(expfac[i], bmfac[i], 1e-4) <<
      "Got unexpected beam factor at pixel number " << i << " for " <<
      pixsize << " arcsec pix band 1";

  //Try with different pixel size
  bm.getBeamFac1(n, 3.0, bmfac);
  const double expfac2[n] = {0.6417129487, 0.8950250709, 1.0, 
			     0.8950250709, 0.6417129487};
  for (unsigned int i = 0; i < n; ++i)
    EXPECT_NEAR(expfac2[i], bmfac[i], 1e-4) <<
      "Got unexpected beam factor at pixel number " << i << 
      " for 3 arcsec pix band 1";

  //Try band 2
  EXPECT_THROW(bm.getBeamFac2(0, pixsize, bmfac), pofdExcept) <<
    "Asking getBeamFac2 for no pixels should throw exception";
  EXPECT_THROW(bm.getBeamFac2(4, pixsize, bmfac), pofdExcept) <<
    "Asking getBeamFac2 for even number of pixels should throw exception";
  EXPECT_THROW(bm.getBeamFac2(n, -pixsize, bmfac), pofdExcept) <<
    "Asking getBeamFac2 for negative pixel size should throw exception";
  EXPECT_THROW(bm.getBeamFac2(n, pixsize, NULL), pofdExcept) <<
    "Null factor array should throw exception";
  bm.getBeamFac2(n, pixsize, bmfac);
  const double expfac_b2[n] = {0.5, 0.8408964152537, 1.0, 0.8408964152537, 0.5};
  for (unsigned int i = 0; i < n; ++i)
    EXPECT_NEAR(expfac_b2[i], bmfac[i], 1e-4) <<
      "Got unexpected beam factor at pixel number " << i << " for " <<
      pixsize << " arcsec pix band 2";

}

TEST(beam2DTest, GetBeam) {
  const double fwhm1 = 15;
  const double fwhm2 = 20;
  const double pixsize = 5.0;
  doublebeam bm(fwhm1, fwhm2);

  const unsigned int n = 5;
  double bmarr[n * n];

  EXPECT_THROW(bm.getBeam1(0, pixsize, bmarr), pofdExcept) <<
    "Asking getBeam1 for no pixels should throw exception";
  EXPECT_THROW(bm.getBeam1(4, pixsize, bmarr), pofdExcept) <<
    "Asking getBeam1 for even number of pixels should throw exception";
  EXPECT_THROW(bm.getBeam1(n, -pixsize, bmarr), pofdExcept) <<
    "Asking getBeam1 for negative pixel size should throw exception";
  EXPECT_THROW(bm.getBeam1(n, pixsize, NULL), pofdExcept) <<
    "Null factor array should throw exception band 1";
  bm.getBeam1(n, pixsize, bmarr);
  const double expfac1[n] = {0.29163225989, 0.734867246, 1.0, 
			     0.734867246, 0.29163225989};
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
    EXPECT_NEAR(expfac1[i] * expfac1[j], bmarr[i * n + j], 1e-4) <<
      "Got unexpected beam value in band 1 at pixel number " << i << " " << j;

  EXPECT_THROW(bm.getBeam2(0, pixsize, bmarr), pofdExcept) <<
    "Asking getBeam2 for no pixels should throw exception";
  EXPECT_THROW(bm.getBeam2(4, pixsize, bmarr), pofdExcept) <<
    "Asking getBeam2 for even number of pixels should throw exception";
  EXPECT_THROW(bm.getBeam2(n, -pixsize, bmarr), pofdExcept) <<
    "Asking getBeam2 for negative pixel size should throw exception";
  EXPECT_THROW(bm.getBeam2(n, pixsize, NULL), pofdExcept) <<
    "Null factor array should throw exception band 1";
  bm.getBeam2(n, pixsize, bmarr);
  const double expfac2[n] = {0.5, 0.8408964152537, 1.0, 0.8408964152537, 0.5};
  for (unsigned int i = 0; i < n; ++i)
    for (unsigned int j = 0; j < n; ++j)
    EXPECT_NEAR(expfac2[i] * expfac2[j], bmarr[i * n + j], 1e-4) <<
      "Got unexpected beam value in band 2 at pixel number " << i << " " << j;

}

////////////////////////////////////
// numberCounts

//Basic Instantiation
TEST(model1DTest, Init) {
  const std::string modelfile("testdata/test1D.txt");
  const unsigned int nknots = 5;
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const double knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

  numberCounts model(modelfile);

  EXPECT_TRUE(model.isValid()) << "Model test case should be valid";
  EXPECT_EQ(5U, model.getNKnots()) << "Unexpected number of model knots";
  EXPECT_NEAR(knotpos[0], model.getMinKnotPosition(), 1e-5) <<
    "Unexpected minimum knot position";
  EXPECT_NEAR(knotpos[nknots-1], model.getMaxKnotPosition(), 1e-5) <<
    "Unexpected maximum knot position";
  for (unsigned int i = 0; i < nknots; ++i)
    EXPECT_NEAR(knotpos[i], model.getKnotPosition(i), 1e-5) <<
      "Unexpected knot position at index " << i;
  for (unsigned int i = 0; i < nknots; ++i)
    EXPECT_NEAR(knotval[i], model.getLog10KnotValue(i), 1e-5) <<
      "Unexpected knot value at index " << i;
  
}

//Number of sources
TEST(model1DTest, Counts) {
  const std::string modelfile("testdata/test1D.txt");
  numberCounts model(modelfile);
  const unsigned int nknots = 5;
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const double knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };

  EXPECT_NEAR(3.06005e6, model.getBaseN0(), 1e4) <<
    "Unexpected base number of sources per area";

  //First check right at knots -- have to skip highest knot where model
  // evaluates to zero.
  for (unsigned int i = 0; i < nknots-1; ++i)
    EXPECT_NEAR(knotval[i], log10(model.getdNdS(knotpos[i])), 1e-3) <<
      "Unexpected counts at knot " << i;

  //Now try some off knot positions
  const unsigned int ntest = 3;
  const double testpos[ntest] = { 0.003, 0.0075, 0.030 };
  const double testval[ntest] = { 8.67247885195, 5.8300749986, 1.07518749639 };
  for (unsigned int i = 0; i < ntest; ++i)
    EXPECT_NEAR(testval[i], log10(model.getdNdS(testpos[i])), 1e-3) <<
      "Unexpected counts at position " << testpos[i];
}

//Flux density per area
TEST(model1DTest, FluxPerArea) {
  const std::string modelfile("testdata/test1D.txt");
  numberCounts model(modelfile);

  EXPECT_NEAR(7233.11, model.getBaseFluxPerArea(), 0.1) <<
    "Unexpected flux per area";
  EXPECT_NEAR(17.7339, model.getBaseFluxSqPerArea(), 0.1) <<
    "Unexpected flux^2 per area";
}

//R testing, single value version.
//These are based on comparisons with IDL code that computes R
TEST(model1DTest, RScalar) {
  const std::string modelfile("testdata/test1D.txt");
  const double fwhm = 15.0;
  const double pixsize = 5.0;
  const double nfwhm = 3.5;
  const unsigned int nbins = 50;

  numberCounts model(modelfile);
  beam bm(fwhm);

  ASSERT_TRUE(model.isValid()) << "Model should be valid";
  EXPECT_THROW(model.getR(0.5, bm, -1.0, nfwhm, nbins), pofdExcept) <<
    "Asking for negative pixel size should throw exception";
  EXPECT_THROW(model.getR(0.5, bm, pixsize, -1.0, nbins), pofdExcept) <<
    "Asking for negative nfwhm should throw exception";
  EXPECT_THROW(model.getR(0.5, bm, pixsize, nfwhm, 0), pofdExcept) <<
    "Asking for no bins should throw exception";

  const unsigned int ntest = 13;
  const double testx[ntest] = {0.001, 0.002, 0.003, 0.004, 0.005,
			       0.006, 0.007, 0.008, 0.009, 0.010,
			       0.015, 0.020, 0.030};
  const double expR[ntest] = {5765.3663, 31000.585, 1467.3891, 175.78777, 
			      35.731806, 10.958402, 4.1296338, 2.0301393,
			      1.1860079, 0.75482242, 0.12561533, 0.019924714,
			      2.2936243e-5};
  double rval, reldiff;
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(testx[i], bm, pixsize, nfwhm, nbins);
    reldiff = fabs((rval - expR[i]) / expR[i]);
    EXPECT_NEAR(0.0, reldiff, 1e-3) << "Unexpectedly large relative R" << 
      "difference for x value " << testx[i] << "; expected " <<
      expR[i] << " got " << rval;
  }

  //Test dependence on nbins -- shouldn't be too sensitive
  const unsigned int nbins2 = nbins * 2;
  for (unsigned int i = 0; i < ntest; ++i) {
    rval = model.getR(testx[i], bm, pixsize, nfwhm, nbins2);
    reldiff = fabs((rval - expR[i]) / expR[i]);
    EXPECT_NEAR(0.0, reldiff, 1e-3) << "Unexpectedly large relative R" << 
      "difference for x value " << testx[i] << "; expected " <<
      expR[i] << " got " << rval << " in nbins doubling test";
  }
}

//R testing, array version
TEST(model1DTest, RArray) {
  const std::string modelfile("testdata/test1D.txt");
  const double fwhm = 15.0;
  const double pixsize = 5.0;
  const double nfwhm = 3.5;
  const unsigned int nbins = 50;

  numberCounts model(modelfile);
  beam bm(fwhm);

  const unsigned int ntest = 10;
  // Computed using IDL code
  const double testx[ntest] = {0.001, 0.002, 0.003, 0.004, 0.005,
			       0.006, 0.007, 0.008, 0.009, 0.010};
  const double expR[ntest] = {5765.3663, 31000.585, 1467.3891, 175.78777, 
			      35.731806, 10.958402, 4.1296338, 2.0301393,
			      1.1860079, 0.75482242};
  double rval[ntest];

  ASSERT_TRUE(model.isValid()) << "Model should be valid";
  EXPECT_THROW(model.getR(ntest, 0.001, 0.010, bm, -1.0, nfwhm, nbins, rval), 
	       pofdExcept) << "Asking for negative pixel size should "
			   << "throw exception";
  EXPECT_THROW(model.getR(ntest, 0.001, 0.010, bm, pixsize, -1.0, nbins, rval), 
	       pofdExcept) << "Asking for negative nfwhm should "
			   << "throw exception";
  EXPECT_THROW(model.getR(ntest, 0.001, 0.010, bm, pixsize, nfwhm, 0, rval), 
	       pofdExcept) << "Asking for no bins should throw exception";

  model.getR(ntest, 0.001, 0.010, bm, pixsize, nfwhm, nbins, rval);
  double reldiff;
  for (unsigned int i = 0; i < ntest; ++i) {
    reldiff = fabs((rval[i] - expR[i]) / expR[i]);
    EXPECT_NEAR(0.0, reldiff, 1e-3) << "Unexpectedly large relative R" << 
      "difference for x value " << testx[i] << "; expected " <<
      expR[i] << " got " << rval[i];
  }
}

////////////////////////////////////
// numberCountsDouble

//Basic Instantiation
TEST(model2DTest, Init) {
  const std::string modelfile("testdata/test2D.txt");
  const unsigned int nknots = 5;
  const double knotpos[nknots] = { 0.002, 0.005, 0.010, 0.020, 0.040 };
  const double knotval[nknots] = { 10.0, 7.0, 5.0, 4.0, -1.0 };
  const unsigned int nsigmas = 3;
  const double sigmapos[nsigmas] = { 0.003, 0.005, 0.020 };
  const double sigmaval[nsigmas] = { 0.1, 0.15, 0.2 };
  const unsigned int noffsets = 1;
  const double offsetpos[noffsets] = { 0.010 };
  const double offsetval[noffsets] = { -0.5 };

  numberCountsDouble model(modelfile);

  EXPECT_TRUE(model.isValid()) << "Model test case should be valid";
  EXPECT_EQ(nknots, model.getNKnots()) << "Unexpected number of model knots";
  EXPECT_EQ(nsigmas, model.getNSigmaKnots()) << 
    "Unexpected number of sigma knots";
  EXPECT_EQ(noffsets, model.getNOffsetKnots()) << 
    "Unexpected number of offset knots";
  EXPECT_EQ(nknots + nsigmas + noffsets, model.getNTotalKnots()) <<
    "Unexpected total number of knots";

  // Band 1 model
  for (unsigned int i = 0; i < nknots; ++i)
    EXPECT_NEAR(knotpos[i], model.getKnotPosition(i), 1e-5) <<
      "Unexpected knot position at index " << i;
  for (unsigned int i = 0; i < nknots; ++i)
    EXPECT_NEAR(knotval[i], model.getLog10KnotValue(i), 1e-5) <<
      "Unexpected knot value at index " << i;
  for (unsigned int i = 0; i < nknots-1; ++i) //Highest knot evaluates to 0
    EXPECT_NEAR(knotval[i], log10(model.getBand1dNdS(knotpos[i])), 1e-3) <<
        "Unexpected band 1 counts from getBand1dNdS at f1 " << knotpos[i];

  // Sigma
  for (unsigned int i = 0; i < nsigmas; ++i)
    EXPECT_NEAR(sigmapos[i], model.getSigmaKnotPosition(i), 1e-5) <<
      "Unexpected sigma position at index " << i;
  for (unsigned int i = 0; i < nsigmas; ++i)
    EXPECT_NEAR(sigmaval[i], model.getSigmaKnotValue(i), 1e-5) <<
      "Unexpected sigma value from getSigmaKnotValue at index " << i;
  for (unsigned int i = 0; i < nsigmas; ++i)
    EXPECT_NEAR(sigmaval[i], model.getSigma(sigmapos[i]), 1e-5) <<
      "Unexpected sigma value from getSigma at f1 " << sigmapos[i];
  
  // Offset
  for (unsigned int i = 0; i < noffsets; ++i)
    EXPECT_NEAR(offsetpos[i], model.getOffsetKnotPosition(i), 1e-5) <<
      "Unexpected offset position at index " << i;
  for (unsigned int i = 0; i < noffsets; ++i)
    EXPECT_NEAR(offsetval[i], model.getOffsetKnotValue(i), 1e-5) <<
      "Unexpected offset value from getOffsetKnotValue at index " << i;
  for (unsigned int i = 0; i < noffsets; ++i)
    EXPECT_NEAR(offsetval[i], model.getOffset(offsetpos[i]), 1e-5) <<
      "Unexpected offset value from getOffsetKnot at f1 " << offsetpos[i];
}

//Number of sources
TEST(model2DTest, Counts) {
  const std::string modelfile("testdata/test2D.txt");

  numberCountsDouble model(modelfile);

  EXPECT_NEAR(3.06005e6, model.getBaseN0(), 1e4) <<
    "Unexpected base number of sources per area";

  //Evaluate at exact value of sigma_i: at f1=0.005 sigma_i is 0.15
  // and offset is -0.5.  The number counts are then the 1D value there
  // (10**7.0) times L(f2 / 0.005; -0.5, 0.15)/0.
  double currf1, sig, off, b1cnts;
  currf1 = 0.005;
  sig = 0.15;
  off = -0.5;
  ASSERT_NEAR(sig, model.getSigma(currf1), 1e-4) <<
    "Unexpected sigma at band 1 flux density " << currf1;
  ASSERT_NEAR(off, model.getOffset(currf1), 1e-4) <<
    "Unexpected offset at band 1 flux density " << currf1;
  b1cnts = 1e7; //Band 1 counts at f1=0.005 
  ASSERT_NEAR(log10(b1cnts), log10(model.getBand1dNdS(currf1)), 1e-3) <<
    "Got unexpected band 1 differential counts at " << currf1;

  // Test for various f2 values
  const unsigned int n2 = 5;
  const double f2[n2] = {0.003, 0.005, 0.006, 0.010, 0.020};
  double expval[n2];
  double frat, expon, cnts;
  for (unsigned int i = 0; i < n2; ++i) {
    frat = f2[i] / currf1;
    expon = (log(frat) - off) / sig;
    cnts = b1cnts * pofd_coverage::isqrt_two_pi / (f2[i] * sig) * 
      exp(-0.5 * expon * expon);
    expval[i] = log10(cnts);
  }
  for (unsigned int i = 0; i < n2; ++i)
    EXPECT_NEAR(expval[i], log10(model.getdNdS(currf1, f2[i])), 1e-3)
  		<< "Unexpected counts for flux 2 value: " << f2[i];
  
  //Now we try a range of f1 values, no longer relying on knowing sig/off
  // ahead of time, but asking the model class for them
  const unsigned int n1 = 5;
  const double f1[n1] = {0.002, 0.003, 0.007, 0.011, 0.03};
  for (unsigned int i = 0; i < n1; ++i) {
    sig = model.getSigma(f1[i]);
    off = model.getOffset(f1[i]);
    b1cnts = model.getBand1dNdS(f1[i]);
    for (unsigned int j = 0; j < n2; ++j) {
      frat = f2[j] / f1[i];
      expon = (log(frat) - off) / sig;
      cnts = b1cnts * pofd_coverage::isqrt_two_pi / (f2[j] * sig) * 
	exp(-0.5 * expon * expon);
      EXPECT_NEAR(log10(cnts), log10(model.getdNdS(f1[i], f2[j])), 1e-3) <<
	"Got unexpected 2 band counts at flux densities " << f1[i] << 
	" " << f2[j];
    }
  }

  //Now compare with an IDL code.  This only supports single sigma/offset
  // values, so we need a different input file
  const std::string modelfile_idl("testdata/test2D_simple.txt");
  numberCountsDouble model_idl(modelfile_idl);

  const unsigned int ntest_idl = 7;
  const double f1_idl[ntest_idl] = {0.002, 0.002, 0.004, 0.006, 0.007,
				    0.010, 0.030};
  const double f2_idl[ntest_idl] = {0.0015, 0.0021, 0.0035, 0.0065, 0.008, 
				    0.011, 0.022};
  const double log10_b1_idl[ntest_idl] = {10, 10, 7.7305874, 6.4739311,
					  6.0291462, 5.0, 1.0751873};
  const double log10_idl[ntest_idl] = {13.126978, 11.158770, 9.9260647,
				       6.8618526, 5.8036657, 5.0157282,
				       3.0796793};
  //First do band 1
  for (unsigned int i = 0; i < ntest_idl; ++i)
    ASSERT_NEAR(log10_b1_idl[i], log10(model_idl.getBand1dNdS(f1_idl[i])), 
		1e-3) << "Got unexpected band 1 differential counts at " 
		      << f1_idl[i];
  
  //Now 2 band
  for (unsigned int i = 0; i < ntest_idl; ++i)
    EXPECT_NEAR(log10_idl[i], log10(model_idl.getdNdS(f1_idl[i], f2_idl[i])),
		1e-3) << "Got unexpected differential counts with "
		      << f1_idl[i] << " " << f2_idl[i];

}

//Flux density per area
TEST(model2DTest, BaseFlux) {
  //Test band 1 fluxes, since they should be the same as
  // for the 1D model in model1DTest::MeanFlux, as test2D
  // is the same model with a colour model appended
  const std::string modelfile1("testdata/test2D.txt");
  numberCountsDouble model1(modelfile1);

  EXPECT_NEAR(7233.11, model1.getBaseFluxPerArea1(), 0.1) <<
    "Unexpected flux per area, band 1";
  EXPECT_NEAR(17.7339, model1.getBaseFluxSqPerArea1(), 0.001) <<
    "Unexpected flux^2 per area, band 1";

  //Now switch over to a simplified model with constant offset and sigma
  const std::string modelfile2("testdata/test2D_simple.txt");
  numberCountsDouble model2(modelfile2);
  double sig = model2.getSigma(0.1);
  ASSERT_NEAR(sig, 0.15, 1e-5) << "Unexpected sigma value";
  double off = model2.getOffset(0.1);
  ASSERT_NEAR(off, -0.4, 1e-5) << "Unexpected offset value";
  double expfac;
  expfac = exp(off + 0.5 * sig * sig);
  EXPECT_NEAR(7233.11 * expfac, model2.getBaseFluxPerArea2(), 0.1) <<
    "Unexpected band 2 flux per area";
  expfac = exp(2 * off + 2.0 * sig * sig);
  EXPECT_NEAR(17.7339 * expfac, model2.getBaseFluxSqPerArea2(), 0.001) <<
    "Unexpected band 2 flux^2 per area";
  
}

//Test R, scalar version, against IDL code
TEST(model2DTest, RScalar) {
  //IDL version supports only constant sigma, offset, so use that
  const std::string modelfile("testdata/test2D_simple.txt");
  numberCountsDouble model(modelfile);
  doublebeam bm(15.0, 20.0);
  const double pixsize = 5.0;
  const double nfwhm = 3.5;
  const unsigned int nbins = 80;

  const unsigned int ntest = 3;
  const double x1[ntest] = {0.001, 0.003, 0.004};
  const double x2[ntest] = {0.0011, 0.004, 0.005};
  const double log10R[ntest] = {7.0781536, 3.5896604, 3.1435221};

  ASSERT_TRUE(model.isValid()) << "Model should be valid";
  EXPECT_THROW(model.getR(x1[0], x2[0], bm, -1.0, nfwhm, nbins), pofdExcept) <<
    "Expected negative pixel size to throw exception";
  EXPECT_THROW(model.getR(x1[0], x2[0], bm, pixsize, -1.0, nbins), 
	       pofdExcept) << "Expected negative nfwhm to throw exception";
  EXPECT_THROW(model.getR(x1[0], x2[0], bm, pixsize, nfwhm, 0), 
	       pofdExcept) << "Expected 0 nbins to throw exception";

  for (unsigned int i = 0; i < ntest; ++i)
    EXPECT_NEAR(log10R[i], log10(model.getR(x1[i], x2[i], bm, pixsize, 
					    nfwhm, nbins)), 1e-3) <<
      "Unexpected R value for x1: " << x1[i] << " x2: " << x2[i];

}


//Test R, array version, against IDL code
TEST(model2DTest, RArray) {
  //IDL version supports only constant sigma, offset, so use that
  const std::string modelfile("testdata/test2D_simple.txt");
  numberCountsDouble model(modelfile);
  doublebeam bm(15.0, 20.0);
  const double pixsize = 5.0;
  const double nfwhm = 3.5;
  const unsigned int nbins = 80;

  const unsigned int n1 = 2;
  const unsigned int n2 = 3;
  const double x1[n1] = {0.001, 0.010};
  const double x2[n2] = {0.0011, 0.004, 0.008};
  const double log10R[n1 * n2] = {  7.078154, 0.94494006, -9.9229534, // x1[0]
				  -28.853694,-0.43926033,  2.3081361}; //x1[1]
  double rarr[n1 * n2];

  //Test failure modes
  ASSERT_TRUE(model.isValid()) << "Model should be valid";
  EXPECT_THROW(model.getR(n1, x1, n2, x2, bm, -1, nfwhm, nbins, rarr), 
	       pofdExcept) << "Expected negative pixel size to throw exception";
  EXPECT_THROW(model.getR(n1, x1, n2, x2, bm, pixsize, -1, nbins, rarr), 
	       pofdExcept) << "Expected negative nfwhm to throw exception";
  EXPECT_THROW(model.getR(n1, x1, n2, x2, bm, pixsize, nfwhm, 0, rarr), 
	       pofdExcept) << "Expected zero bins to throw exception";

  //Test R values
  model.getR(n1, x1, n2, x2, bm, pixsize, nfwhm, nbins, rarr);
  for (unsigned int i = 0; i < n1; ++i)
    for (unsigned int j = 0; j < n2; ++j)
      EXPECT_NEAR(log10R[i * n2 + j], log10(rarr[i * n2 + j]), 1e-3) <<
	"Unexpected R value for x1: " << x1[i] << " x2: " << x2[j];

}


////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running model tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
