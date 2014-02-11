#!/usr/bin/env python3

from __future__ import print_function

import .broken_norm

import math
from collections import OrderedDict

def evaluate_sim1D(filenames, donorm=True):
    """Does the actual evaluation of the likelihood, and returns
    a dictionary of information about the model and results."""
    
    normobj = broken_norm.find_norm1D(filenames, domap=donorm)

    retdict = {}
    retdict['filenames'] = filenames
    retdict['nsims'] = normobj.simset.data.ntotsims
    retdict['nsets'] = normobj.simset.data.nsets
    retdict['n0'] = normobj.simset.data.n0
    retdict['fftsize'] = normobj.simset.data.fftsize
    retdict['fwhm'] = normobj.simset.data.fwhm
    retdict['pixsize'] = normobj.simset.data.pixsize
    retdict['beam_area_pix'] = normobj.simset.data.beam_area_pix
    retdict['beamsq_area_pix'] = normobj.simset.data.beamsq_area_pix
    retdict['area'] = normobj.simset.data.area
    retdict['sigma_inst'] = normobj.simset.data.sigma_inst
    retdict['additional_smoothing'] = normobj.simset.data.additional_smoothing
    if normobj.simset.data.additional_smoothing:
        retdict['extra_smoothing'] = normobj.simset.data.extra_smoothing

    if donorm:
        normvals = normobj.fit()
        retdict['norm'] = normvals[0]
        retdict['norm_plus'] = normvals[1]
        retdict['norm_minus'] = normvals[2]
    
        # Get the variances, etc. as well for the best norm
        stats = normobj.simset.data.get_stats(normvals[0])
        retdict['n0_mean'] = stats[0]
        scat = math.sqrt(stats[1])
        scat_sig = math.sqrt(stats[2])
        retdict['n0_scat'] = scat
        retdict['n0_scat_uncert'] = 0.5 * scat_sig / scat
        retdict['n0_bias'] = normobj.simset.data.n0 - stats[0]
        retdict['n0_like_scat'] = math.sqrt(stats[3])
        retdict['n0_like_scat_uncert'] = 0.5*math.sqrt(stats[4]/stats[3])

    else:
        stats = normobj.simset.data.get_bestn0_stats()
        retdict['n0_mean'] = stats[0]
        scat = math.sqrt(stats[1])
        scat_sig = math.sqrt(stats[2])
        retdict['n0_scat'] = scat
        retdict['n0_scat_uncert'] = 0.5 * scat_sig / scat
        retdict['n0_bias'] = normobj.simset.data.n0 - stats[0]

    return retdict


def evaluate_sim2D(filenames, donorm=True):
    """Does the actual evaluation of the likelihood, and returns
    a dictionary of information about the model and results."""
    
    normobj = broken_norm.find_norm2D(filenames, domap=donorm)

    retdict = {}
    retdict['filenames'] = filenames
    retdict['nsims'] = normobj.simset.data.ntotsims
    retdict['nsets'] = normobj.simset.data.nsets
    retdict['n0'] = normobj.simset.data.n0
    retdict['fftsize'] = normobj.simset.data.fftsize
    retdict['fwhm1'] = normobj.simset.data.fwhm1
    retdict['fwhm2'] = normobj.simset.data.fwhm2
    retdict['pixsize'] = normobj.simset.data.pixsize
    retdict['beam_area_pix1'] = normobj.simset.data.beam_area_pix1
    retdict['beam_area_pix2'] = normobj.simset.data.beam_area_pix2
    retdict['area'] = normobj.simset.data.area
    retdict['sigma_inst1'] = normobj.simset.data.sigma_inst1
    retdict['sigma_inst2'] = normobj.simset.data.sigma_inst2
    retdict['additional_smoothing'] = normobj.simset.data.additional_smoothing
    if normobj.simset.data.additional_smoothing:
        retdict['extra_smoothing1'] = normobj.simset.data.extra_smoothing1
        retdict['extra_smoothing2'] = normobj.simset.data.extra_smoothing2

    if donorm:
        normvals = normobj.fit()
        retdict['norm'] = normvals[0]
        retdict['norm_plus'] = normvals[1]
        retdict['norm_minus'] = normvals[2]
        
        # Get the variances, etc. as well for the best norm
        stats = normobj.simset.data.get_stats(normvals[0])
        retdict['n0_mean'] = stats[0]
        scat = math.sqrt(stats[1])
        scat_sig = math.sqrt(stats[2])
        retdict['n0_scat'] = scat
        retdict['n0_scat_uncert'] = 0.5 * scat_sig / scat
        retdict['n0_bias'] = normobj.simset.data.n0 - stats[0]
        retdict['n0_like_scat'] = math.sqrt(stats[3])
        retdict['n0_like_scat_uncert'] = 0.5*math.sqrt(stats[4]/stats[3])
    else:
        stats = normobj.simset.data.get_bestn0_stats()
        retdict['n0_mean'] = stats[0]
        scat = math.sqrt(stats[1])
        scat_sig = math.sqrt(stats[2])
        retdict['n0_scat'] = scat
        retdict['n0_scat_uncert'] = 0.5 * scat_sig / scat
        retdict['n0_bias'] = normobj.simset.data.n0 - stats[0]

    return retdict

if __name__ == "__main__":
    import argparse
    import pickle
    import glob
    import os.path

    desc = """Fit a set of coverage simulations to determine normalization"""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('filename',action='store',
                        help='Test file listing files for each set')
    parser.add_argument('outfile',action='store',
                        help='Output file')
    parser.add_argument('-d','--datadir',action='store', 
                        help='Directory to look for data files in',
                        default=None)
    msgstr = 'Previous fit results read from here, and new results appended'
    parser.add_argument('-i','--inputfile',action='store',
                        default=None, help=msgstr)
    parser.add_argument('-n', '--nonorm', action='store_true',
                        default=False, help="Don't compute normalization stats")
    parser.add_argument('-t','--twod', action='store_true',
                        default=False, help="Analyze 2D fit")
    parser.add_argument('-v','--verbose',action='store_true', default=False,
                        help='Print informational messages')

    results = parser.parse_args() #Runs on sys.argv by default

    #Input checks
    if not results.datadir is None:
        if not os.path.isdir(results.datadir):
            if not os.path.exists(results.datadir):
                errstr = "datadir %s does not exist"
            else:
                errstr = "datadir %s exists but is not a directory"
            raise IOError(errstr % results.datadir)

    if not os.path.exists(results.filename):
        errstr = "Input file %s not found"
        raise IOError(errstr % results.filename)

    #Possibly read in previous results file
    if not results.inputfile is None:
        if not os.path.exists(results.inputfile):
            errstr = "--inputfile %s not found"
            raise IOError(errstr % results.inputfile)

        with open(results.inputfile, 'rb') as fl:
            fitinfo = pickle.load(fl)

        if not isinstance(fitinfo, OrderedDict):
            fitinfo = OrderedDict(fitinfo)

    else:
        fitinfo = OrderedDict()

    #Read in the file of fits to process
    with open(results.filename, 'r') as fl:
        sets = []
        files = []
        for line in fl.readlines():
            ln = line.strip()
            if len(ln) == 0:
                continue
            if ln[0] == "#":
                continue
            lspl = ln.split()
            sets.append(lspl[0]) #dataset number
            files.append(lspl[1:]) #Files associated, glob patterns

    nsets = len(sets)    
    if nsets == 0:
        raise ValueError("No data sets to process from %s" % results.filename)

    #Main loop
    donorm = not results.nonorm
    for i in range(nsets):
        if results.verbose:
            percent = 100.0 * (i+1) / nsets
            print("Doing set %3i [%5.1f%%]: %s" % (i+1,percent,sets[i]))

            if sets[i] in fitinfo:
                print("Warning: will overwrite set info for %s" % sets[i])

        #Get filenames using glob
        filelist = []
        for fl in files[i]:
            if not results.datadir is None:
                globseq = os.path.join(results.datadir, fl)
            else:
                globseq = fl
            fls = glob.glob(globseq)
            if len(fls) == 0:
                raise ValueError("Found no files to process from %s" % fl)
            filelist.append(fls)

        filelist = [item for sublist in filelist for item in sublist]
        

        if results.twod:
            fitinfo[sets[i]] = evaluate_sim2D(filelist, donorm=donorm)
        else:
            fitinfo[sets[i]] = evaluate_sim1D(filelist, donorm=donorm)

    #Pickle results
    with open(results.outfile, 'wb') as output:
        pickle.dump(fitinfo, output)
        
        
