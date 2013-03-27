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
  for (unsigned int i = 0; i < ntest-1; ++i)
    EXPECT_NEAR(testval[i], log10(model.getdNdS(testpos[i])), 1e-3) <<
      "Unexpected counts at position " << testpos[i];
}

//Flux density per area
TEST(model1DTest, MeanFlux) {
  const std::string modelfile("testdata/test1D.txt");
  numberCounts model(modelfile);

  EXPECT_NEAR(7233.11, model.getMeanFluxPerArea(), 0.1) <<
    "Unexpected mean flux per area";
  EXPECT_NEAR(17.7339, model.getMeanFluxSqPerArea(), 0.1) <<
    "Unexpected mean flux^2 per area";
}

//R testing, single value version
TEST(model1DTest, RSingle) {
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

  const unsigned int ntest = 4;
  const double testx[ntest] = {0.002, 0.003, 0.004, 0.015};
  const double expR[ntest] = {31000.6, 1467.39, 175.788, 0.125615};
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

  const unsigned int ntest = 3;
  const double testx[ntest] = {0.002, 0.003, 0.004};
  double rval[ntest];

  ASSERT_TRUE(model.isValid()) << "Model should be valid";
  EXPECT_THROW(model.getR(ntest, 0.002, 0.004, bm, -1.0, nfwhm, nbins, rval), 
	       pofdExcept) << "Asking for negative pixel size should "
			   << "throw exception";
  EXPECT_THROW(model.getR(ntest, 0.002, 0.004, bm, pixsize, -1.0, nbins, rval), 
	       pofdExcept) << "Asking for negative nfwhm should "
			   << "throw exception";
  EXPECT_THROW(model.getR(ntest, 0.002, 0.004, bm, pixsize, nfwhm, 0, rval), 
	       pofdExcept) << "Asking for no bins should throw exception";

  const double expR[ntest] = {31000.6, 1467.39, 175.788};
  model.getR(ntest, 0.002, 0.004, bm, pixsize, nfwhm, nbins, rval);
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

  for (unsigned int i = 0; i < nknots; ++i)
    EXPECT_NEAR(knotpos[i], model.getKnotPosition(i), 1e-5) <<
      "Unexpected knot position at index " << i;
  for (unsigned int i = 0; i < nknots; ++i)
    EXPECT_NEAR(knotval[i], model.getLog10KnotValue(i), 1e-5) <<
      "Unexpected knot value at index " << i;
  
  for (unsigned int i = 0; i < nsigmas; ++i)
    EXPECT_NEAR(sigmapos[i], model.getSigmaKnotPosition(i), 1e-5) <<
      "Unexpected sigma position at index " << i;
  for (unsigned int i = 0; i < nsigmas; ++i)
    EXPECT_NEAR(sigmaval[i], model.getSigmaKnotValue(i), 1e-5) <<
      "Unexpected sigma value at index " << i;
  
  for (unsigned int i = 0; i < noffsets; ++i)
    EXPECT_NEAR(offsetpos[i], model.getOffsetKnotPosition(i), 1e-5) <<
      "Unexpected offset position at index " << i;
  for (unsigned int i = 0; i < noffsets; ++i)
    EXPECT_NEAR(offsetval[i], model.getOffsetKnotValue(i), 1e-5) <<
      "Unexpected offset value at index " << i;
  

}

//Number of sources
TEST(model2DTest, Counts) {
  const std::string modelfile("testdata/test2D.txt");
  numberCountsDouble model(modelfile);

  EXPECT_NEAR(3.06005e6, model.getBaseN0(), 1e4) <<
    "Unexpected base number of sources per area";
}


////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running model tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
