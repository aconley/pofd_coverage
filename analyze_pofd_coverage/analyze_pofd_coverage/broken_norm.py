from __future__ import print_function

import numpy as np
import scipy.optimize as op
import math

class simdata1D:
    """Class holding the results of a pofd_coverage simulation, 1D case"""

    def __init__(self, filename):
        """Reads in the fits file filename and stores the results"""
        import os.path

        # Figure out how to process based on file extension
        splfl = os.path.splitext(filename)
        if len(splfl) < 2:
            raise IOError("Didn't find extension for {0:s}".format(filename))
        extension = splfl[1][1:].strip().lower()
        if extension == 'fits':
            self._readfits(filename)
        elif extension == 'h5' or extension == 'hdf5':
            self._readhdf5(filename)
        else:
            raise IOError("Didn't recogize extension {0:s}".format(extension))

    def _readfits(self, filename):
        """ Initialize instance from FITS file"""
        import astropy.io.fits
        try:
            hdulist = astropy.io.fits.open(filename)
            if len(hdulist) < 2:
                hdulist.close()
                raise IOError("File doesn't have the expected 2 extensions")
        except OSError as ose:
            errmsg = "Failed to open {0:s} with error: {1:s}"
            raise OSError(errmsg.format(filename, ose.args[0]))

        # This is really just a ton of boring file processing
        try:
            self.n0 = hdulist[1].header['N0']
            self.base_n0 = hdulist[1].header['BASEN0']
            self.fwhm = hdulist[1].header['FWHM']
            self.pixsize = hdulist[1].header['PIXSIZE']
            self.n1 = hdulist[1].header['N1']
            self.n2 = hdulist[1].header['N2']
            self.area = hdulist[1].header['AREA']
            self.sigma_inst = hdulist[1].header['SIGI']
            self.additional_smoothing = hdulist[1].header['ADDSMTH']
            if self.additional_smoothing:
                self.extra_smoothing = hdulist[1].header['ESMOOTH']
                self.sigma_inst_smoothed = hdulist[1].header['SIGISM']
                self.final_fwhm = math.sqrt(self.fwhm**2 + 
                                            self.extra_smoothing**2 )
                self.final_sigma_inst = self.sigma_inst_smoothed
            else:
                self.final_fwhm = self.fwhm
                self.final_sigma_inst = self.sigma_inst
            self.nsims = hdulist[1].header['NSIMS']
            self.fftsize = hdulist[1].header['FFTSIZE']
            self.pofd_delta_version = hdulist[1].header['VERSION'].strip()
            self.date_file_created = hdulist[1].header['DATE']
            if 'NOINIRNG' in hdulist[1].header:
                self.n0_initial_range = hdulist[1].header['N0INIRNG']

            #Beam area
            self.beam_area_pix = hdulist[1].header['BMAREA']
            self.beamsq_area_pix = hdulist[1].header['BMAREASQ']

        except KeyError as ke:
            hdulist.close()
            errmsg = "Error processing header: %s for file %s" % \
                (ke.args,filename)
            raise IOError(errmsg)
        

        #Data
        try :
            if hdulist[1].header['MAPLIKE']:
                self.maplike = True
                self.best_n0 = hdulist[1].data.field('BEST_N0')
                self.best_loglike = hdulist[1].data.field('BEST_LOGLIKE')

                self.min_n0 = hdulist[1].data.field('MIN_N0')
                self.delta_n0 = hdulist[1].data.field('DELTA_N0')
                self.loglike = hdulist[1].data.field('LOGLIKE')
                self.nlike =  hdulist[1].header['NLIKE']
            else:
                self.maplike = False
                self.best_n0 = hdulist[1].data.field('BEST_N0')
                self.best_loglike = hdulist[1].data.field('BEST_LOGLIKE')
        finally:
            hdulist.close()

    def _readhdf5(self, filename):
        """ Initialize instance from HDF5 file"""
        import h5py
        try:
            h = h5py.File(filename, 'r')
        except OSError as ose:
            errmsg = "Failed to open {0:s} with error: {1:s}"
            raise OSError(errmsg.format(filename, ose.args[0]))

        # This is really just a ton of boring ass file processing
        try:
            self.fwhm = h['Beam'].attrs['FWHM'][0]
            self.nfwhm = h['Beam'].attrs['NFWHMKept'][0]
            self.nbeambins = h['Beam'].attrs['NBeamBins'][0]
            self.beam_area_pix = h['Beam'].attrs['BeamPixArea'][0]
            self.beamsq_area_pix = h['Beam'].attrs['BeamSqPixArea'][0]

            self.modeltype = h['Model'].attrs['ModelType'][0].decode
            self.base_n0 = h['Model'].attrs['BaseN0'][0]
            self.model_knotpos = h['Model/KnotPositions'][:]
            self.model_knotvals = h['Model/Log10KnotValues'][:]

            self.n0 = h['Simulations'].attrs['N0'][0]
            self.pixsize = h['Simulations'].attrs['PixelSize'][0]
            self.n1 = h['Simulations'].attrs['N1'][0]
            self.n2 = h['Simulations'].attrs['N2'][0]
            self.area = h['Simulations'].attrs['SimAreaSqDeg'][0]
            self.sigma_inst = h['Simulations'].attrs['SigmaInst'][0]
            self.final_sigma_inst = h['Simulations'].attrs['SigmaInstFinal'][0]

            self.additional_smoothing = bool(h['Simulations'].attrs['Smoothed'][0])
            if self.additional_smoothing:
                self.extra_smoothing = h['Simulations'].attrs['ExtraSmoothing'][0]
            self.nsims = h['Simulations'].attrs['NSims'][0]

            self.fftsize = h['Simulations'].attrs['FFTSize'][0]
            self.n0_initial_range = h['Simulations'].attrs['N0InitRange'][0]

            self.position_type = h['Simulations/PositionGenerator'].attrs['PositionType'][0].decode()
            if self.position_type == 'Clustered':
                self.clustering_k = h['Simulations/PositionGenerator/k'][:]
                self.clustering_pk = h['Simulations/PositionGenerator/Pk'][:]

            self.isHighPass = bool(h['Filter'].attrs['IsHighPassFiltered'][0])
            if self.isHighPass:
                self.highPassFiltScale = h['Filter'].attrs['HighPassFiltScale'][0]
                self.highPassFiltQ = h['Filter'].attrs['HighPassQFactor'][0]

            self.isMatched = bool(h['Filter'].attrs['IsMatchFiltered'][0])
            if self.isMatched:
                self.matchedFWHM = h['Filter'].attrs['MatchedFiltFWHM'][0]
                self.matchedSigmaInst = h['Filter'].attrs['MatchedFiltSigInst'][0]
                self.matchedSigmaConf = h['Filter'].attrs['MatchedFiltSigConf'][0]

        except KeyError as ke:
            h.close()
            errmsg = "Error processing metadata: %s for file %s" % \
                (ke.args, filename)
            raise IOError(errmsg)
        

        #Data
        try :
            self.maplike = bool(h['Simulations'].attrs['DoMapLike'][0])
            self.best_n0 = h['Simulations/BestN0'][:]
            self.best_loglike = h['Simulations/BestLike'][:]
            if self.maplike:
                self.n0likefrac = h['Simulations'].attrs['N0LikeRangeFraction'][0]
                self.min_n0 = h['Simulations/MinN0'][:]
                self.delta_n0 = h['Simulations/DeltaN0'][:]
                self.loglike = h['Simulations/LogLikelihood'][:, :]
        finally:
            h.close()


    def get_likestats(self, normfac):
        """For a given normalization factor (normfac) get statistics
        on the likelihoods.

        Returns a tuple with 2 elements.  The first element is an array
        of the mean probabilities for each simulated image (i.e. an
        estimate of n0), the second is the variance of each simulated image
        (an estimate of the variance in n0).
        """

        if not self.maplike:
            raise LookupError("Likelihood map is not available")
        if normfac <= 0:
            raise ValueError("normfac %f is not valid" % normfac)

        mean = np.zeros(self.nsims)
        var  = np.zeros(self.nsims)
        for like_idx in range(self.nsims):
            curr_like = self.loglike[like_idx,:]
            n0 = self.min_n0[like_idx] + \
                self.delta_n0[like_idx] * np.arange(0, self.nlike)
            maxval = curr_like.max()
            prob = np.exp((curr_like - maxval) / normfac)
            prob *= 1.0 / prob.sum()
            mean[like_idx] = (n0 * prob).sum()
            var[like_idx] = (n0**2 * prob).sum() - mean[like_idx]**2
        
        return (mean,var)
            
    def close_to_edge(self, mindist):
        """Check to see if any likelihood maps are too close to the edge.

        Return True if any in this object are within mindist of the minimum
        or maximum value of n0.  This uses the simpler test of just
        looking at the largest likelihood rather than weighting.
        """
        
        if not self.maplike:
            raise LookupError("Likelihood map is not available")
        if mindist <= 0:
            raise ValueError("Invalid (non-positive) value of mindist: %f" %
                             mindist)

        maxidx = self.loglike.argmax(1)
        mingap = mindist / self.min_n0 #Gap in steps
        if (maxidx < mingap).any() :
            return True
        if (maxidx > (self.nlike - mingap)).any():
            return True
        return False


class simset1D:
    """This reprents a set of simulations with the same parameters."""

    def __init__(self, filenames, delete_edge=True, nsigma=3., domap=True):
        """Reads in the data files in filenames.

       If delete_edge is set, then any input files with entries
        too close to the edge of the likelihood maps are deleted.
        If it is not set, then the code raises an exception if such
        a case is encountered.  nsigma controls how close a fit has
        to be to the edge to be a problem.
        """

        #This is a fairly complex constructor
        if len(filenames) == 0:
            raise ValueError("No filenames provided")
        self.filenames = filenames

        nsig = float(nsigma)
        if nsig <= 0:
            raise ValueError("Invalid nsigma %f" % nsig)

        #The read
        if isinstance(filenames, str):
            self.data = [simdata1D(filenames)]
            self.nsets = 1
        else:
            self.data = [simdata1D(fl) for fl in filenames]
        
        self.nsets = len(self.data)
        self.ntotsims = int(np.array([dat.nsims for dat in self.data]).sum())

        #Make sure they all have the same parameters and store the means
        # as attributes
        val = np.array([dat.n0 for dat in self.data])
        if val.var() / self.data[0].n0 >= 1e-3 :
            raise ValueError("Too much variance in N0 values")
        self.n0 = val.mean()

        val = np.array([dat.fftsize for dat in self.data], dtype=np.float32)
        if val.var() / self.data[0].fftsize >= 1e-3 :
            raise ValueError("Too much variance in FFTSIZE values")
        self.fftsize = int(val.mean())

        val = np.array([dat.fwhm for dat in self.data])
        if val.var() / self.data[0].fwhm >= 1e-3 :
            raise ValueError("Too much variance in FWHM values")
        self.fwhm = val.mean()

        val = np.array([dat.pixsize for dat in self.data])
        if val.var() / self.data[0].pixsize >= 1e-3 :
            raise ValueError("Too much variance in pixel size values")
        self.pixsize = val.mean()

        self.beam_area_pix = self.data[0].beam_area_pix
        self.beamsq_area_pix = self.data[0].beamsq_area_pix

        val = np.array([dat.area for dat in self.data])
        if val.var() / self.data[0].area >= 1e-3 :
            raise ValueError("Too much variance in area values")
        self.area = val.mean()

        val = np.array([dat.sigma_inst for dat in self.data])
        if val.var() / self.data[0].sigma_inst >= 1e-3 :
            raise ValueError("Too much variance in sigma_inst values")
        self.sigma_inst = val.mean()

        val = np.array([dat.final_sigma_inst for dat in self.data])
        if val.var() / self.data[0].final_sigma_inst >= 1e-3 :
            raise ValueError("Too much variance in final_sigma_inst values")
        self.final_sigma_inst = val.mean()

        val = np.array([dat.additional_smoothing for dat in self.data])
        has_smooth = val.any()
        if has_smooth:
            if not val.all():
                errstr = "Some but not all data files had additional smoothing"
                raise ValueError(errstr)
            self.additional_smoothing = True

            val = np.array([dat.extra_smoothing for dat in self.data])
            if val.var() / self.data[0].extra_smoothing >= 1e-3 :
                raise ValueError("Too much variance in extra_smoothing values")
            self.extra_smoothing = val.mean()

        else:
            self.additional_smoothing = False
        
        #See if we have likelihood maps.  If we have any, we require
        # all to have them (or else we eliminate those that don't,
        # depending on if delete_edge is set)
        val = np.array([dat.maplike for dat in self.data])
        if domap and val.any():
            if not val.all():
                wbad = (val == False).nonzero()
                badfiles = ' '.join(np.array(self.filenames)[wbad])
                raise ValueError("Some files don't have likelhood map: %s" %
                                 badfiles)

            #Now we search for files that are too close to the edge.
            # We need an initial guess at the sigma first, so do one with
            # normfac = 1 over all data
            meann0, varn0, varvarn0 = self.get_bestn0_stats()
            stdevn0 = math.sqrt(varn0)
            edge_problem = [dat.close_to_edge(nsig * stdevn0) for dat 
                            in self.data]
            if True in edge_problem:
                if not False in edge_problem:
                    raise ValueError("All likelihood sets had elements near ends")
                if delete_edge:
                    newdat = [tup[0] for tup in zip(self.data, edge_problem)
                              if not tup[1]]
                    msg = "%i simulation sets of %i were eliminated for "+\
                        "being too close to the edge"
                    print(msg % (len(newdat),self.nsets))
                    self.data = newdat
                    self.nsets = len(self.data)
                    self.ntotsims = int(np.array([dat.nsims 
                                                  for dat in self.data]).sum())

                else:
                    badfiles = ' '.join([tup[0] for tup in 
                                         zip(self.filename, edge_problem)
                                         if tup[1]])
                    errstr = "Encountered data sets with edge problems %s"
                    raise ValueError(errstr % badfiles)

            self.maplike = True
        else:
            self.maplike = False



    def get_bestn0_stats(self):
        """Gets statistics of the fits for the best n0

        The return value is a tuple: (mean n0, var n0, var var n0)
        mean n0 is the mean value of the best fit n0 to each P(D)
        var n0 is the variance in n0
        var var n0 is the variance of the estimator for var n0
        """
        
        best_n0 = np.concatenate([dat.best_n0 for dat in self.data])

        #Get 1st, second, 4th moments
        nelem = len(best_n0)
        if nelem < 2:
            raise ValueError("Not enough elements to compute statistics")
        normfac = 1.0 / (nelem - 1.0)

        mom1_n0 = best_n0.mean()
        delta_var = (best_n0 - mom1_n0)**2
        mom2_n0 = normfac * delta_var.sum()
        delta_var **= 2.0
        mom4_n0 = normfac * delta_var.sum()
        var_var_n0 = \
            (mom4_n0 - (nelem-3.0) * mom2_n0**2  / (nelem-1.0))
        var_var_n0 /= nelem

        return (mom1_n0, mom2_n0, var_var_n0)

    def get_stats(self, normfac):
        """Gets statistics of the likelihood maps
        for the specified normalization.

        The return value is a tuple:
        (mean n0, var n0, var var n0, mean V, var V, var var V)
        mean n0 is the mean value of the best fit n0 to each P(D)
        var n0 is the variance in n0
        var var n0 is the variance of the estimator for var n0
        mean V is the mean of the variance estimates for each simulation
        var V is the variance of the variance estimates for each simulation
        var var V is the variance of the variance of the variance estimates (!)

        Each simulation produces two pieces of information: an estimate
        for n0, and an estimate for the variance in n0, V[n0] (for a 
        given normfac).

        The first three elements are based on n0, and the latter three
        on V[n0].
        """
        
        if not self.maplike:
            raise ValueError("Can't use get_stats on data with no likelihood map")

        nfac = float(normfac)
        if nfac <= 0.0:
            raise ValueError("Invalid (non-positive) normalization factor %f" %
                             nfac)

        mn_n0 = []
        var_n0 = []
        for dat in self.data:
            tup = dat.get_likestats(normfac)
            mn_n0.append(tup[0])
            var_n0.append(tup[1])

        #Collapse into single arrays
        mn_n0 = np.concatenate(mn_n0)
        var_n0 = np.concatenate(var_n0)
        
        #Get 1st, second, 4th moments of each
        nelem = len(mn_n0)
        if nelem < 2:
            raise ValueError("Not enough elements in likelihood")
        normfac = 1.0 / (nelem - 1.0)

        mom1_mn_n0 = mn_n0.mean()
        delta_var = (mn_n0 - mom1_mn_n0)**2
        mom2_mn_n0 = normfac * delta_var.sum()
        delta_var **= 2.0
        mom4_mn_n0 = normfac * delta_var.sum()
        var_var_mn_n0 = \
            (mom4_mn_n0 - (nelem-3.0) * mom2_mn_n0**2  / (nelem-1.0))
        var_var_mn_n0 /= nelem

        mom1_var_n0 = var_n0.mean()
        delta_var = (var_n0 - mom1_var_n0)**2
        mom2_var_n0 = normfac * delta_var.sum()
        delta_var **= 2.0
        mom4_var_n0 = normfac * delta_var.sum()
        var_var_var_n0 = \
            (mom4_var_n0 - (nelem-3.0) * mom2_var_n0**2  / (nelem-1.0))
        var_var_var_n0 /= nelem
        
        return (mom1_mn_n0, mom2_mn_n0, var_var_mn_n0,
                mom1_var_n0, mom2_var_n0, var_var_var_n0)
        

class norm_minimizer1D:
    """This is a class for finding the normalization of a bunch of 
    simulations"""

    def __init__(self, filenames, delete_edge=True, nsigma=3., domap=True):
        """Reads in the data files in filenames.

        If delete_edge is set, then any input files with entries
        too close to the edge of the likelihood maps are deleted.
        If it is not set, then the code raises an exception if such
        a case is encountered.  nsigma controls how close a fit has
        to be to the edge to be a problem.
        """
        self.data = simset1D(filenames, delete_edge=delete_edge,
                             nsigma=nsigma, domap=domap)

        #This is the number of sigma away from the mean
        # we try to hit in our variance matching.  Set to 0 to get the
        # best estimate of the normalization, +1 to get the +1 sigma
        # value, etc.
        self._mintype = 0.0

    @property
    def mintype(self):
        return self._mintype

    @mintype.setter
    def mintype(self, value):
        self._mintype = float(value)

    def __call__(self, normfac):
        """Returns the squared difference between the variance in the
        best fit values and the median of the variance estimated for
        each likelihood map.

        This is the objective function to minimize in order to
        determine the normalization.  The only argument is the
        normalization factor."""

        stats = self.data.get_stats(normfac)

        #This is the variance we are trying to match --
        # the actual scatter of the n0 estimates.  However,
        # there is something slightly tricky here in that we
        # allow the user to specify their target with mintype
        if self._mintype == 0:
            targ_var = stats[1]
        else:
            sdev = math.sqrt(stats[2])
            targ_var = stats[1] + self._mintype * sdev
            if targ_var <= 0:
                errstr = "Invalid target variance from: %f + %f * %f "
                raise ValueError(errstr % (stats[1], self._mintype, sdev))

        #This is the variance estimate for this normalization --
        # the mean of the variances estimated for each simulation
        this_var = stats[3]

        return (targ_var - this_var)**2

class find_norm1D:
    """Finds the normalization for the simulations represented by filenames"""
    
    def __init__(self, filenames, nsigma=3.0, domap=True):
        """Loads in data"""
        self.simset = norm_minimizer1D(filenames, True, nsigma, domap=domap)

    def fit(self, verbose=False, get_uncertanties=True):
        if verbose:
            print("Starting normalization fit")
    
        #Get normalization estimate
        norm = op.fminbound(self.simset, 0.1, 2000)

        if get_uncertanties:
            #Do +1 sigma estimate
            if verbose:
                print(" Starting +1 sigma normalization fit")
            self.simset.mintype = 1.0
            norm_plus = op.fminbound(self.simset, 0.1, 2000)

            #And -1 sigma
            if verbose:
                print(" Starting -1 sigma normalization fit")
            self.simset.mintype = -1.0
            norm_minus = op.fminbound(self.simset, 0.1, 2000)
        
            return (norm, norm_plus - norm, norm_minus - norm)
        else:
            return norm

#############################################

class simdata2D:
    """Class holding the results of a simulation, 2D case"""

    def __init__(self, filename):
        """Reads in the fits file filename and stores the results"""
        import os.path

        # Figure out how to process based on file extension
        splfl = os.path.splitext(filename)
        if len(splfl) < 2:
            raise IOError("Didn't find extension for {0:s}".format(filename))
        extension = splfl[1][1:].strip().lower()
        if extension == 'fits':
            self._readfits(filename)
        elif extension == 'h5' or extension == 'hdf5':
            self._readhdf5(filename)
        else:
            raise IOError("Didn't recogize extension {0:s}".format(extension))

    def _readfits(self, filename):
        """ Initialize instance from FITS file"""
        import astropy.io.fits

        try:
            hdulist = astropy.io.fits.open(filename)
            if len(hdulist) < 2:
                hdulist.close()
                raise IOError("File doesn't have the expected 2 extensions")
        except OSError as ose:
            errmsg = "Failed to open {0:s} with error: {1:s}"
            raise OSError(errmsg.format(filename, ose.args[0]))

        #Header params
        try:
            self.n0 = hdulist[1].header['N0']
            self.fwhm1 = hdulist[1].header['FWHM1']
            self.fwhm2 = hdulist[1].header['FWHM2']
            self.pixsize = hdulist[1].header['PIXSIZE']
            self.n1 = hdulist[1].header['N1']
            self.n2 = hdulist[1].header['N2']
            self.area = hdulist[1].header['AREA']
            self.sigma_inst1 = hdulist[1].header['SIGI_1']
            self.sigma_inst2 = hdulist[1].header['SIGI_2']
            self.additional_smoothing = hdulist[1].header['ADDSMTH']
            if self.additional_smoothing:
                self.extra_smoothing1 = hdulist[1].header['ESMOOTH1']
                self.extra_smoothing2 = hdulist[1].header['ESMOOTH2']
                self.sigma_inst_smoothed1 = hdulist[1].header['SIGISM1']
                self.sigma_inst_smoothed2 = hdulist[1].header['SIGISM2']
                self.final_fwhm1 = math.sqrt(self.fwhm1**2 + 
                                             self.extra_smoothing1**2 )
                self.final_fwhm2 = math.sqrt(self.fwhm2**2 + 
                                             self.extra_smoothing2**2 )
                self.final_sigma_inst1 = self.sigma_inst_smoothed1
                self.final_sigma_inst2 = self.sigma_inst_smoothed2
            else:
                self.final_fwhm1 = self.fwhm1
                self.final_fwhm2 = self.fwhm2
                self.final_sigma_inst1 = self.sigma_inst1
                self.final_sigma_inst2 = self.sigma_inst2
            self.nsims = hdulist[1].header['NSIMS']
            self.fftsize = hdulist[1].header['FFTSIZE']
            self.pofd_delta_version = hdulist[1].header['VERSION'].strip()
            self.date_file_created = hdulist[1].header['DATE']
            if 'N0INIRNG' in hdulist[1].header:
                self.n0_initial_range = hdulist[1].header['N0INIRNG']

            #Beam area
            self.beam_area_pix1 = hdulist[1].header['BMAREA1']
            self.beam_area_pix2 = hdulist[1].header['BMAREA2']
            self.beamsq_area_pix1 = hdulist[1].header['BMARESQ1']
            self.beamsq_area_pix2 = hdulist[1].header['BMARESQ2']

        except KeyError as ke:
            hdulist.close()
            errmsg = "Error processing header: %s for file %s" % \
                (ke.args,filename)
            raise IOError(errmsg)
        

        #Data
        try :
            if hdulist[1].header['MAPLIKE']:
                self.maplike = True
                self.best_n0 = hdulist[1].data.field('BEST_N0')
                self.best_loglike = hdulist[1].data.field('BEST_LOGLIKE')

                self.min_n0 = hdulist[1].data.field('MIN_N0')
                self.delta_n0 = hdulist[1].data.field('DELTA_N0')
                self.loglike = hdulist[1].data.field('LOGLIKE')
                self.nlike =  hdulist[1].header['NLIKE']
            else:
                self.maplike = False
                self.best_n0 = hdulist[1].data.field('BEST_N0')
                self.best_loglike = hdulist[1].data.field('BEST_LOGLIKE')
        finally:
            hdulist.close()

    def _readhdf5(self, filename):
        """ Initialize instance from HDF5 file"""
        import h5py
        try:
            h = h5py.File(filename, 'r')
        except OSError as ose:
            errmsg = "Failed to open {0:s} with error: {1:s}"
            raise OSError(errmsg.format(filename, ose.args[0]))

        # This is really just a ton of boring ass file processing
        try:
            self.fwhm1 = h['Beam'].attrs['FWHM1'][0]
            self.fwhm1 = h['Beam'].attrs['FWHM2'][0]
            self.nfwhm = h['Beam'].attrs['NFWHMKept'][0]
            self.nbeambins = h['Beam'].attrs['NBeamBins'][0]
            self.beam_area_pix1 = h['Beam'].attrs['BeamEffArea1'][0]
            self.beam_area_pix2 = h['Beam'].attrs['BeamEffArea2'][0]

            self.modeltype = h['Model'].attrs['ModelType'][0].decode
            self.base_n0 = h['Model'].attrs['BaseN0'][0]
            self.model_knotpos = h['Model/KnotPositions'][:]
            self.model_knotvals = h['Model/Log10KnotValues'][:]
            self.model_sigma_knotpos = h['Model/SigmaKnotPositions'][:]
            self.model_sigma_knotvals = h['Model/SigmaKnotValues'][:]
            self.model_offset_knotpos = h['Model/OffsetKnotPositions'][:]
            self.model_offset_knotvals = h['Model/OffsetKnotValues'][:]

            self.n0 = h['Simulations'].attrs['N0'][0]
            self.pixsize = h['Simulations'].attrs['PixelSize'][0]
            self.n1 = h['Simulations'].attrs['N1'][0]
            self.n2 = h['Simulations'].attrs['N2'][0]
            self.area = h['Simulations'].attrs['SimAreaSqDeg'][0]
            self.sigma_inst1 = h['Simulations'].attrs['SigmaInst1'][0]
            self.final_sigma_inst1 = h['Simulations'].attrs['SigmaInstFinal1'][0]
            self.sigma_inst2 = h['Simulations'].attrs['SigmaInst2'][0]
            self.final_sigma_inst2 = h['Simulations'].attrs['SigmaInstFinal2'][0]

            self.additional_smoothing = bool(h['Simulations'].attrs['Smoothed'][0])
            if self.additional_smoothing:
                self.extra_smoothing1 = h['Simulations'].attrs['ExtraSmoothing1'][0]
                self.extra_smoothing2 = h['Simulations'].attrs['ExtraSmoothing2'][0]
            self.nsims = h['Simulations'].attrs['NSims'][0]

            self.fftsize = h['Simulations'].attrs['FFTSize'][0]
            self.n0_initial_range = h['Simulations'].attrs['N0InitRange'][0]

            self.position_type = h['Simulations/PositionGenerator'].attrs['PositionType'][0].decode()
            if self.position_type == 'Clustered':
                self.clustering_k = h['Simulations/PositionGenerator/k'][:]
                self.clustering_pk = h['Simulations/PositionGenerator/Pk'][:]

            self.isHighPass1 = bool(h['Filter'].attrs['IsHighPassFiltered1'][0])
            if self.isHighPass1:
                self.highPassFiltScale1 = h['Filter'].attrs['HighPassFiltScale1'][0]
                self.highPassFiltQ1 = h['Filter'].attrs['HighPassQFactor1'][0]

            self.isMatched1 = bool(h['Filter'].attrs['IsMatchFiltered1'][0])
            if self.isMatched1:
                self.matchedFWHM1 = h['Filter'].attrs['MatchedFiltFWHM1'][0]
                self.matchedSigmaInst1 = h['Filter'].attrs['MatchedFiltSigInst1'][0]
                self.matchedSigmaConf1 = h['Filter'].attrs['MatchedFiltSigConf1'][0]
            self.isHighPass2 = bool(h['Filter'].attrs['IsHighPassFiltered2'][0])
            if self.isHighPass2:
                self.highPassFiltScale2 = h['Filter'].attrs['HighPassFiltScale2'][0]
                self.highPassFiltQ2 = h['Filter'].attrs['HighPassQFactor2'][0]

            self.isMatched2 = bool(h['Filter'].attrs['IsMatchFiltered2'][0])
            if self.isMatched2:
                self.matchedFWHM2 = h['Filter'].attrs['MatchedFiltFWHM2'][0]
                self.matchedSigmaInst2 = h['Filter'].attrs['MatchedFiltSigInst2'][0]
                self.matchedSigmaConf2 = h['Filter'].attrs['MatchedFiltSigConf2'][0]

        except KeyError as ke:
            h.close()
            errmsg = "Error processing metadata: %s for file %s" % \
                (ke.args, filename)
            raise IOError(errmsg)
        

        #Data
        try :
            self.maplike = bool(h['Simulations'].attrs['DoMapLike'][0])
            self.best_n0 = h['Simulations/BestN0'][:]
            self.best_loglike = h['Simulations/BestLike'][:]
            if self.maplike:
                self.n0likefrac = h['Simulations'].attrs['N0LikeRangeFraction'][0]
                self.min_n0 = h['Simulations/MinN0'][:]
                self.delta_n0 = h['Simulations/DeltaN0'][:]
                self.loglike = h['Simulations/LogLikelihood'][:, :]
        finally:
            h.close()



    def get_likestats(self, normfac):
        """For a given normalization factor (normfac) get statistics
        on the likelihoods.

        Returns a tuple with 2 elements.  The first element is an array
        of the mean probabilities for each simulated image (i.e. an
        estimate of n0), the second is the variance of each simulated image
        (an estimate of the variance in n0).
        """

        if not self.maplike:
            raise LookupError("Likelihood map is not available")
        if normfac <= 0:
            raise ValueError("normfac %f is not valid" % normfac)

        mean = np.zeros(self.nsims)
        var  = np.zeros(self.nsims)
        for like_idx in range(self.nsims):
            curr_like = self.loglike[like_idx,:]
            n0 = self.min_n0[like_idx] + \
                self.delta_n0[like_idx] * np.arange(0, self.nlike)
            maxval = curr_like.max()
            prob = np.exp((curr_like - maxval) / normfac)
            prob *= 1.0 / prob.sum()
            mean[like_idx] = (n0 * prob).sum()
            var[like_idx] = (n0**2 * prob).sum() - mean[like_idx]**2
        
        return (mean,var)
            
    def close_to_edge(self, mindist):
        """Check to see if any likelihood maps are too close to the edge.

        Return True if any in this object are within mindist of the minimum
        or maximum value of n0.  This uses the simpler test of just
        looking at the largest likelihood rather than weighting.
        """
        
        if not self.maplike:
            raise LookupError("Likelihood map is not available")
        if mindist <= 0:
            raise ValueError("Invalid (non-positive) value of mindist: %f" %
                             mindist)

        maxidx = self.loglike.argmax(1)
        mingap = mindist / self.min_n0 #Gap in steps
        if (maxidx < mingap).any() :
            return True
        if (maxidx > (self.nlike - mingap)).any():
            return True
        return False


class simset2D:
    """This reprents a set of simulations with the same parameters."""

    def __init__(self, filenames, delete_edge=True, nsigma=3., domap=True):
        """Reads in the data files in filenames.

       If delete_edge is set, then any input files with entries
        too close to the edge of the likelihood maps are deleted.
        If it is not set, then the code raises an exception if such
        a case is encountered.  nsigma controls how close a fit has
        to be to the edge to be a problem.
        """

        #This is a fairly complex constructor
        if len(filenames) == 0:
            raise ValueError("No filenames provided")
        self.filenames = filenames

        nsig = float(nsigma)
        if nsig <= 0:
            raise ValueError("Invalid nsigma %f" % nsig)

        #The read
        if isinstance(filenames, str):
            self.data = [simdata2D(filenames)]
            self.nsets = 1
        else:
            self.data = [simdata2D(fl) for fl in filenames]
        
        self.nsets = len(self.data)
        self.ntotsims = int(np.array([dat.nsims for dat in self.data]).sum())

        #Make sure they all have the same parameters and store the means
        # as attributes
        val = np.array([dat.n0 for dat in self.data])
        if val.var() / self.data[0].n0 >= 1e-3 :
            raise ValueError("Too much variance in N0 values")
        self.n0 = val.mean()

        val = np.array([dat.fftsize for dat in self.data], dtype=np.float32)
        if val.var() / self.data[0].fftsize >= 1e-3 :
            raise ValueError("Too much variance in FFTSIZE values")
        self.fftsize = int(val.mean())

        val = np.array([dat.fwhm1 for dat in self.data])
        if val.var() / self.data[0].fwhm1 >= 1e-3 :
            raise ValueError("Too much variance in FWHM1 values")
        self.fwhm1 = val.mean()
        val = np.array([dat.fwhm2 for dat in self.data])
        if val.var() / self.data[0].fwhm2 >= 1e-3 :
            raise ValueError("Too much variance in FWHM2 values")
        self.fwhm2 = val.mean()

        val = np.array([dat.pixsize for dat in self.data])
        if val.var() / self.data[0].pixsize >= 1e-3 :
            raise ValueError("Too much variance in pixel size values")
        self.pixsize = val.mean()

        self.beam_area_pix1 = self.data[0].beam_area_pix1
        self.beamsq_area_pix1 = self.data[0].beamsq_area_pix1
        self.beam_area_pix2 = self.data[0].beam_area_pix2
        self.beamsq_area_pix2 = self.data[0].beamsq_area_pix2

        val = np.array([dat.area for dat in self.data])
        if val.var() / self.data[0].area >= 1e-3 :
            raise ValueError("Too much variance in area values")
        self.area = val.mean()

        val = np.array([dat.sigma_inst1 for dat in self.data])
        if val.var() / self.data[0].sigma_inst1 >= 1e-3 :
            raise ValueError("Too much variance in sigma_inst1 values")
        self.sigma_inst1 = val.mean()
        val = np.array([dat.sigma_inst2 for dat in self.data])
        if val.var() / self.data[0].sigma_inst2 >= 1e-3 :
            raise ValueError("Too much variance in sigma_inst2 values")
        self.sigma_inst2 = val.mean()

        val = np.array([dat.final_sigma_inst1 for dat in self.data])
        if val.var() / self.data[0].final_sigma_inst1 >= 1e-3 :
            raise ValueError("Too much variance in final_sigma_inst1 values")
        self.final_sigma_inst1 = val.mean()
        val = np.array([dat.final_sigma_inst2 for dat in self.data])
        if val.var() / self.data[0].final_sigma_inst2 >= 1e-3 :
            raise ValueError("Too much variance in final_sigma_inst2 values")
        self.final_sigma_inst2 = val.mean()

        val = np.array([dat.additional_smoothing for dat in self.data])
        has_smooth = val.any()
        if has_smooth:
            if not val.all():
                errstr = "Some but not all data files had additional smoothing"
                raise ValueError(errstr)
            self.additional_smoothing = True

            val = np.array([dat.extra_smoothing1 for dat in self.data])
            if val.var() / self.data[0].extra_smoothing1 >= 1e-3 :
                raise ValueError("Too much variance in extra_smoothing1 values")
            self.extra_smoothing1 = val.mean()
            val = np.array([dat.extra_smoothing2 for dat in self.data])
            if val.var() / self.data[0].extra_smoothing2 >= 1e-3 :
                raise ValueError("Too much variance in extra_smoothing2 values")
            self.extra_smoothing2 = val.mean()
        else:
            self.additional_smoothing = False
        
        #See if we have likelihood maps.  If we have any, we require
        # all to have them (or else we eliminate those that don't,
        # depending on if delete_edge is set)
        val = np.array([dat.maplike for dat in self.data])
        if domap and val.any():
            if not val.all():
                wbad = (val == False).nonzero()
                badfiles = ' '.join(np.array(self.filenames)[wbad])
                raise ValueError("Some files don't have likelhood map: %s" %
                                 badfiles)

            #Now we search for files that are too close to the edge.
            # We need an initial guess at the sigma first, so do one with
            # normfac = 1 over all data
            meann0, varn0, varvarn0 = self.get_bestn0_stats()
            stdevn0 = math.sqrt(varn0)
            edge_problem = [dat.close_to_edge(nsig * stdevn0) for dat 
                            in self.data]
            if True in edge_problem:
                if not False in edge_problem:
                    raise ValueError("All likelihood sets had elements near ends")
                if delete_edge:
                    newdat = [tup[0] for tup in zip(self.data, edge_problem)
                              if not tup[1]]
                    msg = "%i simulation sets of %i were eliminated for "+\
                        "being too close to the edge"
                    print(msg % (len(newdat),self.nsets))
                    self.data = newdat
                    self.nsets = len(self.data)
                    self.ntotsims = int(np.array([dat.nsims 
                                                  for dat in self.data]).sum())

                else:
                    badfiles = ' '.join([tup[0] for tup in 
                                         zip(self.filename, edge_problem)
                                         if tup[1]])
                    errstr = "Encountered data sets with edge problems %s"
                    raise ValueError(errstr % badfiles)

            self.maplike = True
        else:
            self.maplike = False



    def get_bestn0_stats(self):
        """Gets statistics of the fits for the best n0

        The return value is a tuple: (mean n0, var n0, var var n0)
        mean n0 is the mean value of the best fit n0 to each P(D)
        var n0 is the variance in n0
        var var n0 is the variance of the estimator for var n0
        """
        
        best_n0 = np.concatenate([dat.best_n0 for dat in self.data])

        #Get 1st, second, 4th moments
        nelem = len(best_n0)
        if nelem < 2:
            raise ValueError("Not enough elements to compute statistics")
        normfac = 1.0 / (nelem - 1.0)

        mom1_n0 = best_n0.mean()
        delta_var = (best_n0 - mom1_n0)**2
        mom2_n0 = normfac * delta_var.sum()
        delta_var **= 2.0
        mom4_n0 = normfac * delta_var.sum()
        var_var_n0 = \
            (mom4_n0 - (nelem-3.0) * mom2_n0**2  / (nelem-1.0))
        var_var_n0 /= nelem

        return (mom1_n0, mom2_n0, var_var_n0)

    def get_stats(self, normfac):
        """Gets statistics of the likelihood maps
        for the specified normalization.

        The return value is a tuple:
        (mean n0, var n0, var var n0, mean V, var V, var var V)
        mean n0 is the mean value of the best fit n0 to each P(D)
        var n0 is the variance in n0
        var var n0 is the variance of the estimator for var n0
        mean V is the mean of the variance estimates for each simulation
        var V is the variance of the variance estimates for each simulation
        var var V is the variance of the variance of the variance estimates (!)

        Each simulation produces two pieces of information: an estimate
        for n0, and an estimate for the variance in n0, V[n0] (for a 
        given normfac).

        The first three elements are based on n0, and the latter three
        on V[n0].
        """
        
        if not self.maplike:
            raise ValueError("Can't use get_stats on data with no likelihood map")

        nfac = float(normfac)
        if nfac <= 0.0:
            raise ValueError("Invalid (non-positive) normalization factor %f" %
                             nfac)

        mn_n0 = []
        var_n0 = []
        for dat in self.data:
            tup = dat.get_likestats(normfac)
            mn_n0.append(tup[0])
            var_n0.append(tup[1])

        #Collapse into single arrays
        mn_n0 = np.concatenate(mn_n0)
        var_n0 = np.concatenate(var_n0)
        
        #Get 1st, second, 4th moments of each
        nelem = len(mn_n0)
        if nelem < 2:
            raise ValueError("Not enough elements in likelihood")
        normfac = 1.0 / (nelem - 1.0)

        mom1_mn_n0 = mn_n0.mean()
        delta_var = (mn_n0 - mom1_mn_n0)**2
        mom2_mn_n0 = normfac * delta_var.sum()
        delta_var **= 2.0
        mom4_mn_n0 = normfac * delta_var.sum()
        var_var_mn_n0 = \
            (mom4_mn_n0 - (nelem-3.0) * mom2_mn_n0**2  / (nelem-1.0))
        var_var_mn_n0 /= nelem

        mom1_var_n0 = var_n0.mean()
        delta_var = (var_n0 - mom1_var_n0)**2
        mom2_var_n0 = normfac * delta_var.sum()
        delta_var **= 2.0
        mom4_var_n0 = normfac * delta_var.sum()
        var_var_var_n0 = \
            (mom4_var_n0 - (nelem-3.0) * mom2_var_n0**2  / (nelem-1.0))
        var_var_var_n0 /= nelem
        
        return (mom1_mn_n0, mom2_mn_n0, var_var_mn_n0,
                mom1_var_n0, mom2_var_n0, var_var_var_n0)
        

class norm_minimizer2D:
    """This is a class for finding the normalization of a bunch of 
    simulations"""

    def __init__(self, filenames, delete_edge=True, nsigma=3., domap=True):
        """Reads in the data files in filenames.

        If delete_edge is set, then any input files with entries
        too close to the edge of the likelihood maps are deleted.
        If it is not set, then the code raises an exception if such
        a case is encountered.  nsigma controls how close a fit has
        to be to the edge to be a problem.
        """
        self.data = simset2D(filenames, delete_edge=delete_edge,
                             nsigma=nsigma, domap=domap)

        #This is the number of sigma away from the mean
        # we try to hit in our variance matching.  Set to 0 to get the
        # best estimate of the normalization, +1 to get the +1 sigma
        # value, etc.
        self._mintype = 0.0

    @property
    def mintype(self):
        return self._mintype

    @mintype.setter
    def mintype(self, value):
        self._mintype = float(value)

    def __call__(self, normfac):
        """Returns the squared difference between the variance in the
        best fit values and the median of the variance estimated for
        each likelihood map.

        This is the objective function to minimize in order to
        determine the normalization.  The only argument is the
        normalization factor."""

        stats = self.data.get_stats(normfac)

        #This is the variance we are trying to match --
        # the actual scatter of the n0 estimates.  However,
        # there is something slightly tricky here in that we
        # allow the user to specify their target with mintype
        if self._mintype == 0:
            targ_var = stats[1]
        else:
            sdev = math.sqrt(stats[2])
            targ_var = stats[1] + self._mintype * sdev
            if targ_var <= 0:
                errstr = "Invalid target variance from: %f + %f * %f "
                raise ValueError(errstr % (stats[1], self._mintype, sdev))

        #This is the variance estimate for this normalization --
        # the mean of the variances estimated for each simulation
        this_var = stats[3]

        return (targ_var - this_var)**2

class find_norm2D:
    """Finds the normalization for the simulations represented by filenames"""
    
    def __init__(self, filenames, nsigma=3.0, domap=True):
        """Loads in data"""
        self.simset = norm_minimizer2D(filenames, True, nsigma, domap=domap)

    def fit(self, verbose=False, get_uncertanties=True):
        if verbose:
            print("Starting normalization fit")
    
        #Get normalization estimate
        norm = op.fminbound(self.simset, 0.1, 2000)

        if get_uncertanties:
            #Do +1 sigma estimate
            if verbose:
                print(" Starting +1 sigma normalization fit")
            self.simset.mintype = 1.0
            norm_plus = op.fminbound(self.simset, 0.1, 2000)

            #And -1 sigma
            if verbose:
                print(" Starting -1 sigma normalization fit")
            self.simset.mintype = -1.0
            norm_minus = op.fminbound(self.simset, 0.1, 2000)
        
            return (norm, norm_plus - norm, norm_minus - norm)
        else:
            return norm
