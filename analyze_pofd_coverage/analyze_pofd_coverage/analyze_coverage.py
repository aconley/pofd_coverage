#!python

from __future__ import print_function

import math
from collections import OrderedDict


def evaluate_sim1D(filenames, donorm=True):
    """Does the actual evaluation of the likelihood, and returns
    a dictionary of information about the model and results."""

    from analyze_pofd_coverage import find_norm1D

    normobj = find_norm1D(filenames, domap=donorm)

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

    from analyze_pofd_coverage import find_norm2D

    normobj = find_norm2D(filenames, domap=donorm)

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


def printsum(info, twod=False):
    """ Print summary of analysis"""

    # Figure out longest key name
    maxstrlen = max([len(k) for k in info])+1
    maxstrlen = 30 if maxstrlen > 30 else maxstrlen

    # Figure out if norm has been done
    #  It feels like there should be a better way to get
    #  an arbitrary element of a dictionary, but I can't find one
    key1 = next(iter(info.keys()))
    donorm = 'norm' in info[key1]

    if donorm:
        if not twod:
            fstr = "%-{0:d}s".format(maxstrlen)
            fstr += " %7s %7s %6s %6s %7s %6s %7s %5s"
            print(fstr % ("#name", "sigi", "norm", "norm+", "norm-", "n0",
                          "n0scat", "n0bias", "nsims"))
            fstr = "%-{0:d}s".format(maxstrlen)
            fstr += " %7.5f %7.2f %6.2f %6.2f %7g %6.1f %7.2f %5i"
            for key in info:
                print(fstr % (key, info[key]['sigma_inst'],
                              info[key]['norm'], info[key]['norm_plus'],
                              info[key]['norm_minus'], info[key]['n0'],
                              info[key]['n0_scat'],
                              info[key]['n0_bias'], info[key]['nsims']))
        else:
            fstr = "%-{0:d}s".format(maxstrlen)
            fstr += " %7s %7s %7s %6s %6s %7s %6s %7s %5s"
            print(fstr % ("#name", "sigi1", "sigi2", "norm", "norm+", "norm-",
                          "n0", "n0scat", "n0bias", "nsims"))
            fstr = "%-{0:d}s".format(maxstrlen)
            fstr += " %7.5f %7.5f %7.2f %6.2f %6.2f %7g %6.1f %7.2f %5i"
            for key in info:
                print(fstr % (key, info[key]['sigma_inst1'],
                              info[key]['sigma_inst2'],
                              info[key]['norm'], info[key]['norm_plus'],
                              info[key]['norm_minus'], info[key]['n0'],
                              info[key]['n0_scat'],
                              info[key]['n0_bias'], info[key]['nsims']))
    else:
        if not twod:
            fstr = "%-{0:d}s".format(maxstrlen)
            fstr += " %7s %7s %7s %7s %7s %5s"
            print(fstr % ("#name", "sigi", "n0", "n0scat", "n0bias",
                          "n0bias%", "nsims"))
            fstr = "%-{0:d}s".format(maxstrlen)
            fstr += " %7.5f %7g %6.1f %7.2f %7.4f %5i"
            for key in info:
                print(fstr % (key, info[key]['sigma_inst'],
                              info[key]['n0'], info[key]['n0_scat'],
                              info[key]['n0_bias'],
                              info[key]['n0_bias'] / info[key]['n0'] * 100.0,
                              info[key]['nsims']))
        else:
            fstr = "%-{0:d}s".format(maxstrlen)
            fstr += " %7s %7s %7s %6s %7s %6s %5s"
            print(fstr % ("#name", "sigi1", "sigi2", "n0", "n0scat",
                          "n0bias", "n0bias%", "nsims"))
            fstr = "%-{0:d}s".format(maxstrlen)
            fstr += " %7.5f %7.5f %7g %6.1f %7.2f %7.4f %5i"
            for key in info:
                print(fstr % (key, info[key]['sigma_inst1'],
                              info[key]['sigma_inst2'],
                              info[key]['n0'], info[key]['n0_scat'],
                              info[key]['n0_bias'],
                              info[key]['n0_bias'] / info[key]['n0'] * 100.0,
                              info[key]['nsims']))


def do_analysis(filename, nonorm=False, verbose=False, twod=False,
                datadir=None):
    """ Analyze the new results specified in filename"""

    import glob
    if not os.path.exists(filename):
        errstr = "Filename {0:s} not found".format(filename)
        raise IOError(errstr)

    if results.verbose:
        print("Processing new files specified in {0:s}".format(filename))

    with open(filename, 'r') as fl:
        sets = []
        files = []
        for line in fl.readlines():
            ln = line.strip()
            if len(ln) == 0 or ln[0] == "#":
                continue
            lspl = ln.split()
            sets.append(lspl[0])  # dataset number
            files.append(lspl[1:])  # Files associated, glob patterns

    nsets = len(sets)
    if nsets == 0:
        errstr = "No data sets to process from {0:s}".format(filename)
        raise ValueError(errstr)

    # Figure out how many sets we are doing, optimize string output
    ndigits = math.ceil(math.log10(nsets))
    digfmt = "{{0:{0:d}d}}".format(ndigits)

    #Main analysis loop
    donorm = not nonorm
    percstr = "Doing set " + digfmt + " [{1:5.1f}%]: {2:s}"
    fitinfo = OrderedDict()
    for i in range(nsets):
        if verbose:
            percent = 100.0 * (i + 1) / nsets
            print(percstr.format(i + 1, percent, sets[i]))

        if sets[i] in fitinfo:
            print("Warning: will overwrite set info for {0:s}".format(sets[i]))

        #Get filenames using glob
        filelist = []
        for fl in files[i]:
            if not datadir is None:
                globseq = os.path.join(datadir, fl)
            else:
                globseq = fl
            fls = glob.glob(globseq)
            if len(fls) == 0:
                errstr = "Found no files to process from {0:s}".format(fl)
                raise ValueError(errstr)
            filelist.append(fls)

        filelist = [item for sublist in filelist for item in sublist]

        if twod:
            fitinfo[sets[i]] = evaluate_sim2D(filelist, donorm=donorm)
        else:
            fitinfo[sets[i]] = evaluate_sim1D(filelist, donorm=donorm)

    return fitinfo

if __name__ == "__main__":
    import argparse
    import os.path
    import textwrap
    import sys

    desc = """ Analyze and summarize the results of pofd_coverage runs"""

    epi = textwrap.dedent('''
    This can be run in several modes.  The most basic is to provide an
    input file of pofd_coverage results to analyze.  Alternatively,
    the output from such a run can be read in and summarized.

    The format of --filename (the list of pofd_coverage files to process)
    is a list of name glob pairs (one per line) to be read in an analyzed,
    where name is what the results will be referred to, and glob is a file
    system glob pattern specifying all the files to be combined in the
    analysis of that name.  The results of this analysis are printed to
    the screen, and can be serialized to the specified --outfile as a
    python pickle object.  This mode is specified by providing a --filename
    argument.

    The --inputfile is the results of a previous such run saved with --outfile.
    If --filename is not provided, then this file is simply read in and
    printed.  If --inputfile and --filename are provided, the results of
    processing the data from --filename are appended to those in --inputfile.
    ''')

    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument('-f', '--filename', action='store', default=None,
                        help='File listing files for each set')
    parser.add_argument('-o', '--outfile', action='store', default=None,
                        help='Output file if ')
    parser.add_argument('-d', '--datadir', action='store',
                        help='Directory to look for data files in',
                        default=None)
    msgstr = 'Previous fit results read from here, and new results appended'\
             ' if needed'
    parser.add_argument('-i', '--inputfile', action='store',
                        default=None, help=msgstr)
    parser.add_argument('-n', '--nonorm', action='store_true', default=False,
                        help="Don't compute normalization stats")
    parser.add_argument('--noprint', action='store_true',
                        default=False, help="Don't print summary information")
    parser.add_argument('-t', '--twod', action='store_true',
                        default=False, help="Analyze 2D fit")
    parser.add_argument('-v', '--verbose', action='store_true', default=False,
                        help='Print informational messages')

    results = parser.parse_args()  # Runs on sys.argv by default

    # Quick return if nothing to do...
    if results.filename is None and results.inputfile is None:
        sys.exit(0)

    #Input checks
    if not results.datadir is None and not os.path.isdir(results.datadir):
        if not os.path.exists(results.datadir):
            errstr = "datadir {0:s} does not exist"
        else:
            errstr = "datadir {0:s} exists but is not a directory"
        raise IOError(errstr.format(results.datadir))

    #Possibly read in previous results file
    if not results.inputfile is None:
        import pickle
        if not os.path.exists(results.inputfile):
            errstr = "--inputfile {0:s} not found".format(results.inputfile)
            raise IOError(errstr)
        if results.verbose:
            print("Reading previous results from {0:s}".
                  format(results.inputfile))
        with open(results.inputfile, 'rb') as fl:
            fitinfo = pickle.load(fl)
        if not isinstance(fitinfo, OrderedDict):
            fitinfo = OrderedDict(fitinfo)
    else:
        fitinfo = OrderedDict()

    #Read in the file of outputs to process
    if not results.filename is None:
        newinfo = do_analysis(results.filename, results.nonorm,
                              results.verbose, results.twod,
                              results.datadir)
        fitinfo.update(newinfo)

    # Pickle results
    if not results.outfile is None:
        import pickle
        if results.verbose:
            print("Serializing results to {0:s}".format(results.outfile))
        with open(results.outfile, 'wb') as output:
            pickle.dump(fitinfo, output)

    # Print
    if not results.noprint:
        if results.verbose and not results.filename is None:
            print("#" * 60)
        printsum(fitinfo, twod=results.twod)
