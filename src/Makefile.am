bin_PROGRAMS = pofd_coverage_makeCat pofd_coverage_makeSim \
	pofd_coverage_getPD pofd_coverage_runSim

pofd_coverage_makeCat_SOURCES = global_settings.h ran.h pofdExcept.h \
	numberCounts.h numberCountsDouble.h utility.h beam.h doublebeam.h \
	pofdExcept.cc numberCounts.cc numberCountsDouble.cc utility.cc \
	beam.cc doublebeam.cc \
	pofd_coverage_makeCat.cc

pofd_coverage_makeSim_SOURCES = global_settings.h ran.h pofdExcept.h \
	simImage.h numberCounts.h beam.h utility.h \
	pofdExcept.cc simImage.cc numberCounts.cc utility.cc beam.cc \
	pofd_coverage_makeSim.cc

pofd_coverage_getPD_SOURCES = global_settings.h ran.h pofdExcept.h \
	numberCounts.h beam.h PDFactory.h PD.h utility.h \
	pofdExcept.cc numberCounts.cc beam.cc PDFactory.cc PD.cc utility.cc \
	pofd_coverage_getPD.cc

pofd_coverage_runSim_SOURCES = global_settings.h ran.h pofdExcept.h \
	numberCounts.h beam.h PDFactory.h PD.h utility.h simImage.h \
	simManager.h \
	pofdExcept.cc numberCounts.cc beam.cc PDFactory.cc PD.cc utility.cc \
	simImage.cc simManager.cc \
	pofd_coverage_runSim.cc
