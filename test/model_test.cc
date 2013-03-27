//Model testing
#include<iostream>

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

////////////////////////////////////////////

GTEST_API_ int main(int argc, char **argv) {
  std::cout << "Running model tests\n";

  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
