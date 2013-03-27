noinst_LTLIBRARIES = libcover.la

libcover_la_SOURCES = pofdExcept.cc utility.cc numberCounts.cc \
	numberCountsDouble.cc beam.cc doublebeam.cc \
	PDFactory.cc PDFactoryDouble.cc PD.cc PDDouble.cc \
	simImage.cc simImageDouble.cc simManager.cc simManagerDouble.cc

bin_PROGRAMS = pofd_coverage_makeCat pofd_coverage_makeSim \
	pofd_coverage_getPD pofd_coverage_runSim

pofd_coverage_makeCat_SOURCES = pofd_coverage_makeCat.cc
pofd_coverage_makeCat_LDADD = libcover.la

pofd_coverage_makeSim_SOURCES = pofd_coverage_makeSim.cc
pofd_coverage_makeSim_LDADD = libcover.la

pofd_coverage_getPD_SOURCES = pofd_coverage_getPD.cc
pofd_coverage_getPD_LDADD = libcover.la

pofd_coverage_runSim_SOURCES = pofd_coverage_runSim.cc
pofd_coverage_runSim_LDADD = libcover.la
