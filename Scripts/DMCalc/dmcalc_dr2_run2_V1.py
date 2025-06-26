#!/usr/bin/python3
################### WIDEBAND DM ESTIMATION SCRIPT ###################
# Script for estimating the DM of an observation using PSRCHIVE and 
# TEMPO2. This code is extracted from another one written for LOFAR 
# data analysis by Caterina (trimtim.py). 
#
# Dependencies: 
# PSRCHIVE python interface: 
# http://psrchive.sourceforge.net/manuals/python/
# SKLEARN: https://scikit-learn.org/stable/install.html
# TEMPO2: https://bitbucket.org/psrsoft/tempo2
# SCIPY
#
# Usage: 
# ./dmcalc.py -E test.par -M test.sm test.fits
#
# For more options and information, please check help section.
#
# If you have a directory with all the model files in one directory 
# (PWD/templates/) with '.sm' extension and all the parameter files 
# in 'ephemerides' directory with the name 'JXXXX-YYYY.par', you do 
# not need to give the -E and -M options and simply do
#
# ./dmcalc.py test.fits
#
#####################################################################

# import modules...
import subprocess
import os
import sys
import numpy as np
import psrchive
import argparse
import time
import warnings
warnings.filterwarnings("ignore")
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.patches import Rectangle


start = time.time()


def existing_file(filename):
	if not(os.path.exists(filename) and os.path.isfile(filename)):
		raise argparse.ArgumentTypeError("**"+filename+" does not exist**")
	else:
		return filename

parser = argparse.ArgumentParser(description='Code for measuring in-band '+ 
                                 'DM for pulsar data in psrfits format.')
parser.add_argument('files', nargs='+', type=existing_file, 
					help='The list of fits file(s) for processing')
parser.add_argument('-E', '--ephem', type=existing_file, 
					help='Ephemeris file to update the model. Exits if not '+
					      'given or is not available in "PWD/ephemerides" '+
					      'directory')
parser.add_argument('-M', '--model', nargs='+', type=existing_file,
					help='Model template for ToA generation. Exits if not '+ 
					     'given or is not available in "PWD/templates" '+
					     'directory')
parser.add_argument('-nch','--nchan', type=int, default=8,
					help='Number of frequency channels to use while '+
						 'estimating DM (Def: 8)')
parser.add_argument('-b3n','--b3nchan', type=int, default=16,
					help='Number of frequency channels to use in '+ 
					     'band3 uGMRT data (Def: 16)')
parser.add_argument('-b4n','--b4nchan', type=int, default=16,
					help='Number of frequency channels to use in '+ 
					     'band4 uGMRT data (Def: 16)')                         
parser.add_argument('-b5n','--b5nchan', type=int, default=8,
					help='Number of frequency channels to use in '+ 
					     'band5 uGMRT data (Def: 8)')
parser.add_argument('-q', '--quiet', action='store_true', 
							help='Only print warnings')
parser.add_argument('-v5','--v5orold', action='store_true',
					help='Employs backend delay correction for old uGMRT processed '+
						 'data using the measured delays. Def: False')

"""
Estimates the Dispersion Measure (DM) from the data in psrfits file format.

Returns the value of DM with its uncertainty and reduced chi-square from
either tempo2 or an MCMC method based on fitting a second order polynomial
using scipy.optimize.

Parameters
----------
file(s) :  Input file(s) in psrfits format

tempo2  :  bool, optional. default: False. If True, performs DM calculation 
           using  tempo2. Otherwise, uses  the MCMC  method to  estimate DM.
          
ephem   :  Ephemeris (or parameter) file  of the  pulsar. This is  required 
           to update the model. It can be  given as a command line argument. 
           If it is available in "PWD/ephemerides" folder, one can use that.
           Giving the file with this option overrides the default one.

model   :  Template profile for cross-correlating  with the observation  to
           obtain DM. It can be given as a command line argument, otherwise
           it will look  for a matching one in  "PWD/ephemerides" directory
           and if found, will use that instead. One can use this  option to
           override the default selection.
           
nchan   : int, optional, default: 8. Number of frequency channels to use in
          the estimation of DM.

b3nchan : int, optional, default: 16. Number of frequency channels in band3
          of uGMRT data.
          
b5nchan : int, optional, default: 8. Number of frequency channels in band5
          of uGMRT data.

v5       : This is used for processing the uGMRT data folded with pinta v5 
           or earlier versions. The backend delays are measured and are kept
           based on the NCRA internal report.

quiet    : bool, optional,  default: False. Supresses all print  statements
           except warnings and errors.

Returns
-------
Dispersion Measure with uncertainty.

Notes
-----

Examples
--------
# (a) for a simple DM estimation with built-in fitting function and default
# directories:
#
dmcalc.py inputfile.fits
#
# (b) for using tempo2 to measure DM:
#
dmcalc.py inputfile.fits
#
# (c) to use different ephemeris and template files:
#
dmcalc.py -E aaaa.eph -M bbbb.fits inputfile.fits
#

"""
# Module that connects all the others and does the job.
def main():
	
	# parses the input arguments
	args = parser.parse_args()

	# checks status of quiet and tempo2
	quiet=False
	if args.quiet:
		quiet=True
	if args.ephem != None:
		ephemeris = args.ephem
	else:
		ephemeris = "ephemerides/"+ar_psr+".par"
		if not (os.path.exists(ephemeris)):
			print("\nError! No parameter file is given. Exiting!\n")
			sys.exit(1)
	if not quiet:
		print ("Pulsar Parameter file is:"+ephemeris+'\n')
		
	if not quiet:
		print("Loading the archive files for DM estimation... "),
	
	

	# loads the data and template file(s)
	pwd=os.getcwd()
	archives = []
	finalarchives = []
	model = []
	for filename in args.files:
		archives.append(psrchive.Archive_load(filename))
		finalarchives.append(psrchive.Archive_load(filename))
	narch = len(archives)
	for filename in args.model:
		model.append(psrchive.Archive_load(filename))

	print("Please wait...")	

	nmod = len(model)	
	arfrq = np.zeros(narch)
	modfrq = np.zeros(nmod)
	ar_psr = archives[0].get_source()
	# Reads the DM from ephemeris and applies it to both data and template files.
	# This step is no longer needed, but is kept as is for being safe than sorry.
	templ_dm = 0.0
	with open (args.ephem, 'r') as read_eph:
		for line in read_eph:
			if line.startswith('DM\t') or line.startswith('DM '):
				templ_dm = float (line.split()[1])
	for i in range (narch):
		arfrq[i] = archives[i].get_centre_frequency()
		archives[i].set_dispersion_measure(templ_dm)
		archives[i].set_ephemeris(ephemeris)
		archives[i].update_model()
		finalarchives[i].set_ephemeris(ephemeris)
		finalarchives[i].update_model()
	for j in range(nmod):
		modfrq[j] = model[j].get_centre_frequency()
		#model[j].set_dispersion_measure(templ_dm)
		#model[j].set_ephemeris(ephemeris)
		model[j].update_model()

	# If there are more than one input and template file, then the following part
	# makes the ordering of those files correct, if it was in a different order.
	# Additionally, for uGMRT data prior to pinta v6.2, it does the proper time 
	# corrections. It will also make the nchan of the data and template as per the
	# input. 
	if (narch == nmod and narch > 1):
		modpos = np.zeros(nmod)
		for i in range(narch):
			if (300. < archives[i].get_centre_frequency() < 500.):
				if args.v5orold:
					archives[i] = Correct_delay(archives[i])
				archives[i].fscrunch_to_nchan(args.b3nchan)
			if (1100. < archives[i].get_centre_frequency() < 1500.):				
				if args.v5orold:
					archives[i] = Correct_delay(archives[i])
				archives[i].fscrunch_to_nchan(args.b5nchan)
			if (525. < archives[i].get_centre_frequency() < 1000.):
				if args.v5orold:
					archives[i] = Correct_delay(archives[i])
				archives[i].fscrunch_to_nchan(args.b4nchan)
			#if (525. < archives[i].get_centre_frequency() < 1000.):
				 #if args.v5orold:
					#archives[i] = Correct_delay(archives[i])
				#archives[i].fscrunch_to_nchan(args.b4nchan)    
		for i in range(nmod):
			if (300. < model[i].get_centre_frequency() < 500.):
				model[i].fscrunch_to_nchan(args.b3nchan)
			if (1100. < model[i].get_centre_frequency() < 1500.):
				model[i].fscrunch_to_nchan(args.b5nchan)
			if (525. < model[i].get_centre_frequency() < 1000.):
				model[i].fscrunch_to_nchan(args.b4nchan)
           	#	if (525. < model[i].get_centre_frequency() < 1000.):
		#		model[i].fscrunch_to_nchan(args.b4nchan)    
		for i in range(narch):
			for j in range(nmod):
				if (np.around(arfrq[i]) == np.around(modfrq[j])):
					modpos[i] = j

	if not quiet:
		print(" done!")
	if (nmod == 2):
		model[0], model[1] = model[int(modpos[0])], model[int(modpos[1])]
	if (nmod == 3):
		model[0], model[1], model[2] = model[int(modpos[0])], model[int(modpos[1])], model[int(modpos[2])]
	nchan = args.nchan
	if (narch == 1):
		if (args.b3nchan != 16):
			nchan = args.b3nchan
		if (args.b5nchan != 8):
			nchan = args.b5nchan
		if (args.b4nchan != 16):
			nchan = args.b4nchan   
		archives[0].fscrunch_to_nchan(nchan)
		model[0].fscrunch_to_nchan(nchan)


	for i in range(narch):
		get_finalTOA(archives[i], model[i], ephemeris, templ_dm, quiet)



	if not quiet:
		print("Calculating DM.......")

	dm, dmerror, chisq, freq_array_orig, resid_array_orig, residerr_array_orig, filtered_array_fil = dm_estimate(ar_psr,ephemeris)

	if not quiet:
		print("Successfully DM Calculated.")

	
	os.remove(ar_psr+"_allToAs_copied.tim")
	#print(dminfof)
	#print(freq_array_in)
	#print(residual_array_in)
	#print(filtered_array_fil)
	#print(dm)
	print("Please wait...")
	filename_fil = filtered_array_fil[0]
	freqs_fil = filtered_array_fil[1].astype(float)
	toas_fil = filtered_array_fil[2].astype(float)
	toaE_fil = filtered_array_fil[3].astype(float)
	tel_fil = filtered_array_fil[4]

	tempfile = ar_psr+"filtered_toas.txt"
	f = open(tempfile,"w+")
	head="FORMAT 1\n"
	f.write('%s' % head)
	for i in range(np.size(freqs_fil)):
		f.write('%s %.8f %.18f %.6f %s\n' % (filename_fil[i], freqs_fil[i], toas_fil[i], toaE_fil[i], tel_fil[i]))
	f.close()

	#temp1 = os.popen("tempo2 -output general2 -f %s %s -s \"1111111 {freq} {pre} {err}\n\" | grep '1111111'" 
	#				% (ephemeris,tempfile)).read()

	#lines2 = temp1.splitlines()
	#prefit_resid = ar_psr+"prefit_resid.txt" 
	#with open(prefit_resid, "w") as f1:
	#	f1.write('\n'.join(lines2))
	
	
	tempo2_output = os.popen("tempo2 -f %s %s -nofit -fit dm -output general2 -s \"11111 {freq} {pre} {post} {res} {err}e-6\n\" | grep '11111'" 
					% (ephemeris,tempfile)).read()
	lines1 = tempo2_output.splitlines()
	postfit_resid = ar_psr+"_resid.txt" 
	with open(postfit_resid, "w") as f:
		f.write('\n'.join(lines1))
	
	#allTOAtim = ar_psr+"_allToAs.tim"



	#unfil_resid = np.loadtxt(ar_psr+"prefit_unfil_resid.txt")
	prepost_resid = np.loadtxt(ar_psr+"_resid.txt")

	os.remove(tempfile)
	#unfiltered
	orig_freqs = freq_array_orig
	orig_resid = resid_array_orig
	orig_toasE = residerr_array_orig
	

	#filtered(prefit)
	init_freqs = prepost_resid[:, 1]
	init_resid = prepost_resid[:, 2]
	init_toasE = prepost_resid[:, 4]

	#filtered(postfit)
	final_freqs = prepost_resid[:, 1]
	final_resid = prepost_resid[:, 3]
	final_toasE = prepost_resid[:, 4]

	os.remove(ar_psr+"_resid.txt")
	
	init_resid -= np.median(init_resid)
	final_resid -= np.median(final_resid)


	orig_resid =	orig_resid*1e+6 
	init_resid =	init_resid*1e+6  
	final_resid = final_resid*1e+6
	orig_toasE = orig_toasE*1e+6
	init_toasE = init_toasE * 1e+6
	final_toasE = final_toasE*1e+6
	

	prefit_rms = np.zeros(np.size(dm)); postfit_rms = np.zeros(np.size(dm))
	med_toaE = np.zeros(np.size(dm)); centre_freqs = np.zeros(np.size(dm))
	bw = np.zeros(np.size(dm))
	prefit_rms[0] = np.sqrt(np.cov(init_resid, aweights=init_toasE))
	postfit_rms[0] = np.sqrt(np.cov(final_resid, aweights=final_toasE))
	med_toaE[0] = np.median(final_toasE)

	#print(prefit_rms[0])

	
	mjd_start=archives[0].get_Integration(0).get_start_time().in_days()
	mjd_end=archives[0].get_Integration(0).get_end_time().in_days()
	ar_mjd = mjd_start + (mjd_end-mjd_start)/2.
	ar_psr = archives[0].get_source()
	ar_tel = archives[0].get_telescope()
	# Removing the DM and DMEPOCH from a copy of the ephemeris file given.
	oldpar = open(ephemeris,"r")
	partmp = ar_psr+'_tmp.par'
	newpar = open(partmp,"w+")
	for i, line in enumerate(oldpar):
		if not line.lstrip().startswith('DM'):
				if not line.lstrip().startswith('DMEPOCH'):
					newpar.write(line)
	oldpar.close()
	newpar.close()

	# updating the ephemeris file with measured DM
	dmline = "DM             "+str(dm[0])+"\t1\t"+str(dmerror[0])
	dmepochline  = "DMEPOCH	       "+str(round(ar_mjd,2))
	f = open(partmp,'a')
	f.write('%s\n%s\n' % (dmline, dmepochline))
	f.close()
	# Correcting the observed files with the obtained DM and getting the 
	# DM corrected ToAs for plotting and statistics.


	for i in range(narch):
		#archives[i].set_ephemeris(partmp)
		#archives[i].set_dispersion_measure(dm[0])
		#archives[i].update_model()
		finalarchives[i].set_ephemeris(partmp)
		finalarchives[i].set_dispersion_measure(dm[0])
		finalarchives[i].update_model()
		finalarchives[i].tscrunch()
		finalarchives[i].dedisperse()
		finalarchives[i].remove_baseline()


	


	#print(bw)

	# Setting the correct band flags and obtaining the BW, centre freq of all the 
	# files for writing them out to a file.
	if not quiet:
		print("Generating DM time series.......")

	bandflag = [None] * np.size(dm)
	if (len(archives) == 1):
		centre_freqs[0] = archives[0].get_centre_frequency()
		bw[0] = archives[0].get_bandwidth()
		if (300. < centre_freqs[0] < 500.0):
			bandflag[0] = 'band3'
		if (1160. < centre_freqs[0] < 1460.0):
			bandflag[0] = 'band5'
		if (525. < centre_freqs[0] < 1000.0):
			bandflag[0] = 'band4'
	
	if (len(archives) == 2):
		cf1 = archives[0].get_centre_frequency()
		cf2 = archives[1].get_centre_frequency()
		cfs = cf1 + cf2
		if (1500. < cfs < 1900.): 
			bandflag[0] = "band3+5"
			for i in range(len(archives)):
				cfrq = archives[i].get_centre_frequency()
				if (300. < cfrq < 500.):
					bw[1] = archives[i].get_bandwidth()
					centre_freqs[1] = archives[i].get_centre_frequency()
					bandflag[1] = "band3"
					condition = (init_freqs < 500.) & (init_freqs > 300.)
					tx = np.extract(condition,init_resid)
					tz = np.extract(condition,init_toasE)
					prefit_rms[1] = np.sqrt(np.cov(tx, aweights=tz))
					condition = (final_freqs < 500.) & (final_freqs > 300.)
					ty = np.extract(condition,final_resid)
					tz = np.extract(condition,final_toasE)
					postfit_rms[1] = np.sqrt(np.cov(ty, aweights=tz))
					med_toaE[1] = np.median(tz)
				if (1160. < cfrq < 1460.):
					bw[2] = archives[i].get_bandwidth()
					centre_freqs[2] = archives[i].get_centre_frequency()
					bandflag[2] = "band5"
					condition = (init_freqs < 1460.) & (init_freqs > 1160.)
					tx = np.extract(condition,init_resid)
					tz = np.extract(condition,init_toasE)
					prefit_rms[2] = np.sqrt(np.cov(tx, aweights=tz))
					condition = (final_freqs < 1460.) & (final_freqs > 1160.)
					ty2 = np.extract(condition,final_resid)
					tz2 = np.extract(condition,final_toasE)
					postfit_rms[2] = np.sqrt(np.cov(ty2, aweights=tz2))
					med_toaE[2] = np.median(tz2)
				bw[0] = bw[1] + bw[2]		
				centre_freqs[0] = (centre_freqs[1] + centre_freqs[2])/2.
		if (1000. < cfs < 1300.): 
			bandflag[0] = "band3+4"
			for i in range(len(archives)):
				cfrq = archives[i].get_centre_frequency()
				if (300. < cfrq < 500.):
					bw[1] = archives[i].get_bandwidth()
					centre_freqs[1] = archives[i].get_centre_frequency()
					bandflag[1] = "band3"
					condition = (init_freqs < 500.) & (init_freqs > 300.)
					tx = np.extract(condition,init_resid)
					tz = np.extract(condition,init_toasE)
					prefit_rms[1] = np.sqrt(np.cov(tx, aweights=tz))
					condition = (final_freqs < 500.) & (final_freqs > 300.)
					ty = np.extract(condition,final_resid)
					tz = np.extract(condition,final_toasE)
					postfit_rms[1] = np.sqrt(np.cov(ty, aweights=tz))
					med_toaE[1] = np.median(tz)
				if (525. < cfrq < 1000.):
					bw[2] = archives[i].get_bandwidth()
					centre_freqs[2] = archives[i].get_centre_frequency()
					bandflag[2] = "band4"
					condition = (init_freqs < 1000.) & (init_freqs > 525.)
					tx = np.extract(condition,init_resid)
					tz = np.extract(condition,init_toasE)
					prefit_rms[2] = np.sqrt(np.cov(tx, aweights=tz))
					condition = (final_freqs < 1000.) & (final_freqs > 525.)
					ty2 = np.extract(condition,final_resid)
					tz2 = np.extract(condition,final_toasE)
					postfit_rms[2] = np.sqrt(np.cov(ty2, aweights=tz2))
					med_toaE[2] = np.median(tz2)
				bw[0] = bw[1] + bw[2]		
				centre_freqs[0] = (centre_freqs[1] + centre_freqs[2])/2.
		if (2000. < cfs < 2300.): 
			bandflag[0] = "band4+5"
			for i in range(len(archives)):
				cfrq = archives[i].get_centre_frequency()
				if (525. < cfrq < 1000.):
					bw[1] = archives[i].get_bandwidth()
					centre_freqs[1] = archives[i].get_centre_frequency()
					bandflag[1] = "band4"
					condition = (init_freqs < 1000.) & (init_freqs > 525.)
					tx = np.extract(condition,init_resid)
					tz = np.extract(condition,init_toasE)
					prefit_rms[1] = np.sqrt(np.cov(tx, aweights=tz))
					condition = (final_freqs < 1000.) & (final_freqs > 525.)
					ty = np.extract(condition,final_resid)
					tz = np.extract(condition,final_toasE)
					postfit_rms[1] = np.sqrt(np.cov(ty, aweights=tz))
					med_toaE[1] = np.median(tz)
				if (1160. < cfrq < 1460.):
					bw[2] = archives[i].get_bandwidth()
					centre_freqs[2] = archives[i].get_centre_frequency()
					bandflag[2] = "band5"
					condition = (init_freqs < 1460.) & (init_freqs > 1160.)
					tx = np.extract(condition,init_resid)
					tz = np.extract(condition,init_toasE)
					prefit_rms[2] = np.sqrt(np.cov(tx, aweights=tz))
					condition = (final_freqs < 1460.) & (final_freqs > 1160.)
					ty2 = np.extract(condition,final_resid)
					tz2 = np.extract(condition,final_toasE)
					postfit_rms[2] = np.sqrt(np.cov(ty2, aweights=tz2))
					med_toaE[2] = np.median(tz2)
				bw[0] = bw[1] + bw[2]		
				centre_freqs[0] = (centre_freqs[1] + centre_freqs[2])/2.
	
	if (len(archives) > 2):
		bandflag[0] = "band3+4+5"
		for i in range(len(archives)):
			cfrq = archives[i].get_centre_frequency()
			if (300. < cfrq < 500.):
				bw[4] = archives[i].get_bandwidth()
				centre_freqs[4] = archives[i].get_centre_frequency()
				bandflag[4] = "band3"
				condition = (init_freqs < 500.) & (init_freqs > 300.)
				tx = np.extract(condition,init_resid)
				tz = np.extract(condition,init_toasE)
				prefit_rms[4] = np.sqrt(np.cov(tx, aweights=tz))
				condition = (final_freqs < 500.) & (final_freqs > 300.)
				ty = np.extract(condition,final_resid)
				tz = np.extract(condition,final_toasE)
				postfit_rms[4] = np.sqrt(np.cov(ty, aweights=tz))
				med_toaE[4] = np.median(tz)
			if (525. < cfrq < 1000.):
				bw[5] = archives[i].get_bandwidth()
				centre_freqs[5] = archives[i].get_centre_frequency()
				bandflag[5] = "band4"
				condition = (init_freqs < 1000.) & (init_freqs > 525.)
				tx = np.extract(condition,init_resid)
				tz = np.extract(condition,init_toasE)
				prefit_rms[5] = np.sqrt(np.cov(tx, aweights=tz))
				condition = (final_freqs < 1000.) & (final_freqs > 525.)
				ty2 = np.extract(condition,final_resid)
				tz2 = np.extract(condition,final_toasE)
				postfit_rms[5] = np.sqrt(np.cov(ty2, aweights=tz2))
				med_toaE[5] = np.median(tz2)
			if (1160. < cfrq < 1460.):
				bw[6] = archives[i].get_bandwidth()
				centre_freqs[6] = archives[i].get_centre_frequency()
				bandflag[6] = "band5"
				condition = (init_freqs < 1460.) & (init_freqs > 1160.)
				tx = np.extract(condition,init_resid)
				tz = np.extract(condition,init_toasE)
				prefit_rms[6] = np.sqrt(np.cov(tx, aweights=tz))
				condition = (final_freqs < 1460.) & (final_freqs > 1160.)
				ty2 = np.extract(condition,final_resid)
				tz2 = np.extract(condition,final_toasE)
				postfit_rms[6] = np.sqrt(np.cov(ty2, aweights=tz2))
				med_toaE[6] = np.median(tz2)

			for j in range(i+1,len(archives)):
				cfrq1 = archives[j].get_centre_frequency()
				cfs = cfrq + cfrq1
				if (1000. < cfs < 1300.):
					bw[1] = (archives[i].get_bandwidth()+archives[j].get_bandwidth())
					centre_freqs[1] = (archives[i].get_centre_frequency()+archives[j].get_centre_frequency())/2
					bandflag[1] = "band3+4"
					condition = (init_freqs < 1000.) & (init_freqs > 300.)
					tx = np.extract(condition,init_resid)
					tz = np.extract(condition,init_toasE)
					prefit_rms[1] = np.sqrt(np.cov(tx, aweights=tz))
					condition2 = (final_freqs < 1000.) & (final_freqs > 300.)
					ty2 = np.extract(condition2,final_resid)
					tz2 = np.extract(condition2,final_toasE)
					postfit_rms[1] = np.sqrt(np.cov(ty2, aweights=tz2))
					med_toaE[1] = np.median(tz2)

				if (1500. < cfs < 1900.):
					bw[2] = (archives[i].get_bandwidth()+archives[j].get_bandwidth())
					centre_freqs[2] = (archives[i].get_centre_frequency()+archives[j].get_centre_frequency())/2
					bandflag[2] = "band3+5"
					condition = condition = ((1160 < init_freqs) & (init_freqs < 1460.)) | ((300. < init_freqs) & (init_freqs < 500.))
					tx = np.extract(condition,init_resid)
					tz = np.extract(condition,init_toasE)
					prefit_rms[2] = np.sqrt(np.cov(tx, aweights=tz))
					condition2 = ((1160 < final_freqs) & (final_freqs < 1460.)) | ((300. < final_freqs) & (final_freqs < 500.))
					ty2 = np.extract(condition2,final_resid)
					tz2 = np.extract(condition2,final_toasE)
					postfit_rms[2] = np.sqrt(np.cov(ty2, aweights=tz2))
					med_toaE[2] = np.median(tz2)


				if (2000. < cfs < 2300.):
					bw[3] = (archives[i].get_bandwidth()+archives[j].get_bandwidth())
					centre_freqs[3] = (archives[i].get_centre_frequency()+archives[j].get_centre_frequency())/2
					bandflag[3] = "band4+5"
					condition = (init_freqs < 1460.) & (init_freqs > 525.)
					tx = np.extract(condition,init_resid)
					tz = np.extract(condition,init_toasE)
					prefit_rms[3] = np.sqrt(np.cov(tx, aweights=tz))
					condition2 = (final_freqs < 1460.) & (final_freqs > 525.)
					ty2 = np.extract(condition2,final_resid)
					tz2 = np.extract(condition2,final_toasE)
					postfit_rms[3] = np.sqrt(np.cov(ty2, aweights=tz2))
					med_toaE[3] = np.median(tz2)					 	
								




			bw[0] = bw[4] + bw[5] +bw[6]
			centre_freqs[0] = (centre_freqs[4] + centre_freqs[5] + centre_freqs[6])/3.





	# Printing the results to the file and also in the terminal
	f= open(ar_psr+"_DM_timeseries.txt",'a')
	for i in range(np.size(dm)):
		f.write('%.6f %.6f %.6f %.2f %.7f %.7f %.7f %.2f %.2f %s %s\n' %(
			ar_mjd, dm[i], dmerror[i], chisq[i], prefit_rms[i], postfit_rms[i], med_toaE[i], centre_freqs[i], bw[i], ar_tel, bandflag[i]))
	f.close()

	if not quiet:
		print("Successfully DM time series file created.")
	# Creating a par file with DMMODEL parameters
	#dmmodelpar = ar_psr+"_"+str("%.f" % bw[0])+"MHz.DMMODEL.par"
	if not quiet:
		print("Generating DMMODEL and DMX par files.")
	dmmodelpar = ar_psr+".DMMODEL.par"
	if (os.path.isfile(dmmodelpar)):
		oldpar = open(dmmodelpar,"r")
	if not (os.path.isfile(dmmodelpar)):
		f1 = open(ephemeris,"r")
		partmp = dmmodelpar
		dmmodelpar = open(partmp,"w+")
		for line in f1:
			dmmodelpar.write(line)
		dmmodelpar.write("DMMODEL DM 0\n")
		dmmodelpar.close()
		del dmmodelpar
		f1.close()
		#dmmodelpar = ar_psr+"_"+str("%.f" % bw[0])+"MHz.DMMODEL.par"
		dmmodelpar = ar_psr+".DMMODEL.par"
		oldpar = open(dmmodelpar,"r")
	#oldpar = dmmodelpar
	partmp = ar_psr+'_tmp.par'
	newpar = open(partmp,"w+")
	for i, line in enumerate(oldpar):
		dmomjdstr="DMOFF "+str("%.6f" % ar_mjd)
		#print(dmomjdstr)
		if line.lstrip().startswith(dmomjdstr):
			print("\nERROR: DM for MJD %.6f already exists. Maybe try moving the DMMODEL and DMX files and retry." % ar_mjd)
			sys.exit(1)
		if not line.lstrip().startswith('CONSTRAIN'):
					newpar.write(line)
	oldpar.close()
	newpar.close()
	
	ttmp = os.popen("mv %s %s" % (partmp,dmmodelpar)).read()
	f1 = open(dmmodelpar,"a")
	f1.write("DMOFF %.6f %.8f %.8f\n" %(ar_mjd,dm[0]-templ_dm,dmerror[0]))
	f1.write("CONSTRAIN DMMODEL")
	f1.close()

	# Creating a par file with DMX parameters
	#dmxpar = ar_psr+"_"+str("%.f" % bw[0])+".DMX.par"
	dmxpar = ar_psr+".DMX.par"
	if (os.path.isfile(dmxpar)):
		with open(dmxpar,'r') as f:
			last_line = f.readlines()[-1]
			last_dmx =int(last_line.strip().split()[0].split('_')[1])
			this_dmx = last_dmx+1
			dmx1 = f"DMX_{this_dmx:04d}"; dmx2 = f"DMXEP_{this_dmx:04d}"
			dmx3 = f"DMXR1_{this_dmx:04d}"; dmx4 = f"DMXR2_{this_dmx:04d}"
			dmx5 = f"DMXF1_{this_dmx:04d}"; dmx6 = f"DMXF2_{this_dmx:04d}"
		if (len(dm)==7):
			f1 = open(dmxpar,"a")
			f1.write("%s\t%.8e\t0\t%.8e\n" % (dmx1, dm[0]-templ_dm, dmerror[0]))
			f1.write("%s\t%.5f\n" % (dmx2, ar_mjd))
			f1.write("%s\t%.5f\n%s\t%.5f\n" % (dmx3, ar_mjd-0.5, dmx4, ar_mjd+0.5))
			f1.write("%s\t%.5f\n%s\t%.5f\n" % (dmx5, centre_freqs[4], dmx6, centre_freqs[6]))
		if (len(dm)==3):
			f1 = open(dmxpar,"a")
			f1.write("%s\t%.8e\t0\t%.8e\n" % (dmx1, dm[0]-templ_dm, dmerror[0]))
			f1.write("%s\t%.5f\n" % (dmx2, ar_mjd))
			f1.write("%s\t%.5f\n%s\t%.5f\n" % (dmx3, ar_mjd-0.5, dmx4, ar_mjd+0.5))
			f1.write("%s\t%.5f\n%s\t%.5f\n" % (dmx5, centre_freqs[1], dmx6, centre_freqs[2]))
			f1.close()
		if (len(dm)==1):
			f1 = open(dmxpar,"a")
			f1.write("%s\t%.8e\t0\t%.8e\n" % (dmx1, dm[0]-templ_dm, dmerror[0]))
			f1.write("%s\t%.5f\n" % (dmx2, ar_mjd))
			f1.write("%s\t%.5f\n%s\t%.5f\n" % (dmx3, ar_mjd-0.5, dmx4, ar_mjd+0.5))
			f1.write("%s\t%.5f\n%s\t%.5f\n" % (dmx5, centre_freqs[0]-bw[0]/2., dmx6, centre_freqs[0]+bw[0]/2))
			f1.close()
	
		#dmxpar = ar_psr+"_"+str("%.f" % bw[0])+".DMX.par"
	dmxpar = ar_psr+".DMX.par"
	if not (os.path.isfile(dmxpar)):

		if (len(dm) == 7):
			f1 = open(ephemeris,"r")
			partmp = dmxpar
			dmxpar = open(partmp,"w+")
			for line in f1:
				dmxpar.write(line)
				#dmxpar.write("DMX\t\t6.500000\n")
			dmxpar.write("DMX_0001\t%.8e\t0\t%.8e\n" % (dm[0]-templ_dm,dmerror[0]))
			dmxpar.write("DMXEP_0001\t%.5f\n" % (ar_mjd))
			dmxpar.write("DMXR1_0001\t%.5f\n" % (ar_mjd-0.5))
			dmxpar.write("DMXR2_0001\t%.5f\n" % (ar_mjd+0.5))
			dmxpar.write("DMXF1_0001\t%.5f\n" % (centre_freqs[4]))
			dmxpar.write("DMXF2_0001\t%.5f\n" % (centre_freqs[6]))
			dmxpar.close()
			f1.close()

		if (len(dm) == 3):
			f1 = open(ephemeris,"r")
			partmp = dmxpar
			dmxpar = open(partmp,"w+")
			for line in f1:
				dmxpar.write(line)
				#dmxpar.write("DMX\t\t6.500000\n")
			dmxpar.write("DMX_0001\t%.8e\t0\t%.8e\n" % (dm[0]-templ_dm,dmerror[0]))
			dmxpar.write("DMXEP_0001\t%.5f\n" % (ar_mjd))
			dmxpar.write("DMXR1_0001\t%.5f\n" % (ar_mjd-0.5))
			dmxpar.write("DMXR2_0001\t%.5f\n" % (ar_mjd+0.5))
			dmxpar.write("DMXF1_0001\t%.5f\n" % (centre_freqs[1]))
			dmxpar.write("DMXF2_0001\t%.5f\n" % (centre_freqs[2]))
			dmxpar.close()
			f1.close()
		if (len(dm) ==1):
			f1 = open(ephemeris,"r")
			partmp = dmxpar
			dmxpar = open(partmp,"w+")
			for line in f1:
				dmxpar.write(line)
				#dmxpar.write("DMX\t\t6.500000\n")
			dmxpar.write("DMX_0001\t%.8e\t0\t%.8e\n" % (dm[0]-templ_dm,dmerror[0]))
			dmxpar.write("DMXEP_0001\t%.5f\n" % (ar_mjd))
			dmxpar.write("DMXR1_0001\t%.5f\n" % (ar_mjd-0.5))
			dmxpar.write("DMXR2_0001\t%.5f\n" % (ar_mjd+0.5))
			dmxpar.write("DMXF1_0001\t%.5f\n" % (centre_freqs[0]-bw[0]/2))
			dmxpar.write("DMXF2_0001\t%.5f\n" % (centre_freqs[0]+bw[0]/2))
			dmxpar.close()
			f1.close()
	if not quiet:
		print("Parameter files created.")



	if not quiet:
		print("Generating plots....")

	if (len(finalarchives) == 2):
		cf1 = archives[0].get_centre_frequency()
		cf2 = archives[1].get_centre_frequency()
		cfs = cf1 + cf2
		if (1500. < cfs < 1900.): 

			prof2Db3 = []
			profb3 = []
			prof2Db5 = []
			profb5 = []
			b3_bw = []; b5_bw = []; b3_freq = []; b5_freq = []; b3_nbin = []; b5_nbin = []
			for i in range(len(finalarchives)):
				if (finalarchives[i].get_centre_frequency() < 500.):
					ar_nchan = finalarchives[i].get_nchan()
					b3_nbin  = finalarchives[i].get_nbin()
					b3_bw = finalarchives[i].get_bandwidth()
					b3_freq = finalarchives[i].get_centre_frequency()
					prof2Db3 = finalarchives[i].get_data()[:,0,:,:].flatten().reshape(ar_nchan,b3_nbin)
					prof = finalarchives[i].clone()
					prof.fscrunch()
					profb3 = prof.get_data().flatten()
					profb3 /= np.max(profb3)
				if (finalarchives[i].get_centre_frequency() > 1000.):
					ar_nchan = finalarchives[i].get_nchan()
					b5_nbin  = finalarchives[i].get_nbin()
					b5_bw = finalarchives[i].get_bandwidth()
					b5_freq = finalarchives[i].get_centre_frequency()
					prof2Db5 = finalarchives[i].get_data()[:,0,:,:].flatten().reshape(ar_nchan,b5_nbin)
					prof = finalarchives[i].clone()
					prof.fscrunch()
					profb5 = prof.get_data().flatten()
					profb5 /= np.max(profb5)
		
			condition = orig_freqs < 500.
			orig_b3fr = np.extract(condition,orig_freqs)
			orig_b3re = np.extract(condition,orig_resid)
			orig_b3Er = np.extract(condition,orig_toasE)

			condition = orig_freqs > 1000.
			orig_b5fr = np.extract(condition,orig_freqs)
			orig_b5re = np.extract(condition,orig_resid)
			orig_b5Er = np.extract(condition,orig_toasE)

			condition = init_freqs < 500.
			init_b3fr = np.extract(condition,init_freqs)
			init_b3re = np.extract(condition,init_resid)
			init_b3Er = np.extract(condition,init_toasE)

			condition = init_freqs > 1000.
			init_b5fr = np.extract(condition,init_freqs)
			init_b5re = np.extract(condition,init_resid)
			init_b5Er = np.extract(condition,init_toasE)

			condition = final_freqs < 500.
			final_b3fr = np.extract(condition,final_freqs)
			final_b3re = np.extract(condition,final_resid)
			final_b3Er = np.extract(condition,final_toasE)

			condition = final_freqs > 1000.
			final_b5fr = np.extract(condition,final_freqs)
			final_b5re = np.extract(condition,final_resid)
			final_b5Er = np.extract(condition,final_toasE)

			fig = plt.figure(3, figsize=(8, 6))
			fig.subplots_adjust(hspace=0.05)
			ax0 = plt.subplot2grid((9, 8), (0,0), rowspan=3, colspan=3)
			ax1 = plt.subplot2grid((9, 8), (3,0), rowspan=1, colspan=3)
			ax2 = plt.subplot2grid((9, 8), (5,0), rowspan=3, colspan=3)
			ax3 = plt.subplot2grid((9, 8), (8,0), rowspan=1, colspan=3)
		
			ax4 = plt.subplot2grid((9, 8), (0,4), colspan=4, rowspan=3)
			ax5 = plt.subplot2grid((9, 8), (3,4), colspan=4, rowspan=3)
			ax6 = plt.subplot2grid((9, 8), (6,4), colspan=4, rowspan=3)

			#ax7 = plt.subplot2grid((9, 8), (6,4), colspan=4, rowspan=3)

			leg1 = Rectangle((0, 0), 0, 0, alpha=0.0)
		
			ax0.imshow((np.sqrt(prof2Db5**2))**0.5, origin='lower', extent=(0,b5_nbin-1,(np.around(b5_freq)-b5_bw/2),(np.around(b5_freq)+b5_bw/2)), aspect='auto', cmap='hot')
			ax0.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=8)
			ax0.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
			ax1.plot(np.arange(b5_nbin, dtype=float),profb5, color='black', linewidth=0.5)
			ax1.set_xlim(0,b5_nbin-1)
			ax1.set_ylabel('Intensity', fontweight='bold', fontsize=8)

			ax2.imshow((np.sqrt(prof2Db3**2))**0.5, origin='lower', extent=(0,b3_nbin-1,(np.around(b3_freq)-b3_bw/2),(np.around(b3_freq)+b3_bw/2)), aspect='auto', cmap='hot')
			ax2.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=8)
			ax2.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
			ax3.plot(np.arange(b3_nbin, dtype=float),profb3, color='black', linewidth=0.5)
			ax3.set_xlim(0,b3_nbin-1)
			ax3.set_xlabel('Pulse Phase (bins)', fontweight='bold', fontsize=8)
			ax3.set_ylabel('Intensity', fontweight='bold', fontsize=8)

			ax4.errorbar(orig_b3fr, orig_b3re, yerr=orig_b3Er, fmt='.', color='#D81B60', capsize=2)
			ax4.errorbar(orig_b5fr, orig_b5re, yerr=orig_b5Er, fmt='.', color='#1E88E5', capsize=2)
			ax4.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b5_freq)+b3_bw/0.8))
			ax4.grid()
			ax4.legend([leg1], ['Prefit: Unfiltered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
			ax4.axes.xaxis.set_ticklabels([])

			ax5.errorbar(init_b3fr, init_b3re, yerr=init_b3Er, fmt='.', color='#D81B60', capsize=2)
			ax5.errorbar(init_b5fr, init_b5re, yerr=init_b5Er, fmt='.', color='#1E88E5', capsize=2)
			ax5.grid()
			ax5.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b5_freq)+b3_bw/0.8))
			ax5.legend([leg1], ['Prefit: Filtered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
			ax5.axes.xaxis.set_ticklabels([])
			ax5.set_ylabel(r'ToA Residuals ($\mu$s)', fontweight='bold', fontsize=8)
		
			ax6.errorbar(final_b3fr, final_b3re, yerr=final_b3Er, fmt='.', color='#D81B60', capsize=2)
			ax6.errorbar(final_b5fr, final_b5re, yerr=final_b5Er, fmt='.', color='#1E88E5', capsize=2)
			ax6.grid()
			ax6.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b5_freq)+b3_bw/0.8))
			ax6.legend([leg1], ['Postfit'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
			ax6.set_xlabel('Frequency (MHz)', fontweight='bold', fontsize=9)
			
			#ax7.errorbar(x1, y1, yerr=e1, fmt='.', color='#D81B60', capsize=2)
			#ax7.errorbar(x2, y2, yerr=e2, fmt='.', color='green', capsize=2)
			#ax7.errorbar(x3, y3, yerr=e3, fmt='.', color='#1E88E5', capsize=2)
			#ax7.grid()
			#ax7.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b5_freq)+b3_bw/0.8))
			#ax7.legend([leg1], ['Postfit'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
			#ax7.set_xlabel('Frequency (MHz)', fontweight='bold', fontsize=9)				
			#ax7.set_ylabel('P.fit.res($\mu$s))', fontweight='bold', fontsize=9)

			fig.suptitle('Source: PSR %s; MJD: %.4f; Prefit Wrms: %.2f $\mu$s; Postfit Wrms: %.2f $\mu$s\nMedian ToA Err: %.2f $\mu$s; DM: %.6f $\pm$ %.6f pc cm$^{-3}$;  Reduced $\chi^2$: %.2f'%(ar_psr, ar_mjd, prefit_rms[0], postfit_rms[0], med_toaE[0], dm[0], dmerror[0], chisq[0]), fontsize=10, fontweight='bold')


			dirplot=os.path.join(pwd,ar_psr+"_"+ar_tel+"_plots")
			if not os.path.exists(dirplot):
		   		os.makedirs(dirplot)
			plotfile=dirplot+"/"+ar_psr+"_"+str(ar_mjd)+"_"+str(centre_freqs[0])+"_"+ar_tel+"_DMfitResid.pdf"
			plt.savefig(plotfile, format='pdf')
			plt.close()
		if (1000. < cfs < 1300.): 

			prof2Db3 = []
			profb3 = []
			prof2Db4 = []
			profb4 = []
			b3_bw = []; b4_bw = []; b3_freq = []; b4_freq = []; b3_nbin = []; b4_nbin = []
			for i in range(len(finalarchives)):
				if (finalarchives[i].get_centre_frequency() < 500.):
					ar_nchan = finalarchives[i].get_nchan()
					b3_nbin  = finalarchives[i].get_nbin()
					b3_bw = finalarchives[i].get_bandwidth()
					b3_freq = finalarchives[i].get_centre_frequency()
					prof2Db3 = finalarchives[i].get_data()[:,0,:,:].flatten().reshape(ar_nchan,b3_nbin)
					prof = finalarchives[i].clone()
					prof.fscrunch()
					profb3 = prof.get_data().flatten()
					profb3 /= np.max(profb3)
				if (525. < finalarchives[i].get_centre_frequency() < 1000.):
					ar_nchan = finalarchives[i].get_nchan()
					b4_nbin  = finalarchives[i].get_nbin()
					b4_bw = finalarchives[i].get_bandwidth()
					b4_freq = finalarchives[i].get_centre_frequency()
					prof2Db4 = finalarchives[i].get_data()[:,0,:,:].flatten().reshape(ar_nchan,b4_nbin)
					prof = finalarchives[i].clone()
					prof.fscrunch()
					profb4 = prof.get_data().flatten()
					profb4 /= np.max(profb4)
		
			condition = orig_freqs < 500.
			orig_b3fr = np.extract(condition,orig_freqs)
			orig_b3re = np.extract(condition,orig_resid)
			orig_b3Er = np.extract(condition,orig_toasE)

			condition = (orig_freqs < 1000.) & (orig_freqs > 525.)
			orig_b4fr = np.extract(condition,orig_freqs)
			orig_b4re = np.extract(condition,orig_resid)
			orig_b4Er = np.extract(condition,orig_toasE)

			condition = init_freqs < 500.
			init_b3fr = np.extract(condition,init_freqs)
			init_b3re = np.extract(condition,init_resid)
			init_b3Er = np.extract(condition,init_toasE)

			condition = (init_freqs < 1000.) & (init_freqs > 525.)
			init_b4fr = np.extract(condition,init_freqs)
			init_b4re = np.extract(condition,init_resid)
			init_b4Er = np.extract(condition,init_toasE)

			condition = final_freqs < 500.
			final_b3fr = np.extract(condition,final_freqs)
			final_b3re = np.extract(condition,final_resid)
			final_b3Er = np.extract(condition,final_toasE)

			condition = (final_freqs < 1000.) & (final_freqs > 525.)
			final_b4fr = np.extract(condition,final_freqs)
			final_b4re = np.extract(condition,final_resid)
			final_b4Er = np.extract(condition,final_toasE)

			fig = plt.figure(3, figsize=(8, 6))
			fig.subplots_adjust(hspace=0.05)
			ax0 = plt.subplot2grid((9, 8), (0,0), rowspan=3, colspan=3)
			ax1 = plt.subplot2grid((9, 8), (3,0), rowspan=1, colspan=3)
			ax2 = plt.subplot2grid((9, 8), (5,0), rowspan=3, colspan=3)
			ax3 = plt.subplot2grid((9, 8), (8,0), rowspan=1, colspan=3)
		
			ax4 = plt.subplot2grid((9, 8), (0,4), colspan=4, rowspan=3)
			ax5 = plt.subplot2grid((9, 8), (3,4), colspan=4, rowspan=3)
			ax6 = plt.subplot2grid((9, 8), (6,4), colspan=4, rowspan=3)

			#ax7 = plt.subplot2grid((9, 8), (6,4), colspan=4, rowspan=3)

			leg1 = Rectangle((0, 0), 0, 0, alpha=0.0)
		
			ax0.imshow((np.sqrt(prof2Db4**2))**0.5, origin='lower', extent=(0,b4_nbin-1,(np.around(b4_freq)-b4_bw/2),(np.around(b4_freq)+b4_bw/2)), aspect='auto', cmap='hot')
			ax0.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=8)
			ax0.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
			ax1.plot(np.arange(b4_nbin, dtype=float),profb4, color='black', linewidth=0.5)
			ax1.set_xlim(0,b4_nbin-1)
			ax1.set_ylabel('Intensity', fontweight='bold', fontsize=8)

			ax2.imshow((np.sqrt(prof2Db3**2))**0.5, origin='lower', extent=(0,b3_nbin-1,(np.around(b3_freq)-b3_bw/2),(np.around(b3_freq)+b3_bw/2)), aspect='auto', cmap='hot')
			ax2.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=8)
			ax2.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
			ax3.plot(np.arange(b3_nbin, dtype=float),profb3, color='black', linewidth=0.5)
			ax3.set_xlim(0,b3_nbin-1)
			ax3.set_xlabel('Pulse Phase (bins)', fontweight='bold', fontsize=8)
			ax3.set_ylabel('Intensity', fontweight='bold', fontsize=8)

			ax4.errorbar(orig_b3fr, orig_b3re, yerr=orig_b3Er, fmt='.', color='#D81B60', capsize=2)
			ax4.errorbar(orig_b4fr, orig_b4re, yerr=orig_b4Er, fmt='.', color='green', capsize=2)
			ax4.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b4_freq)+b3_bw/0.8))
			ax4.grid()
			ax4.legend([leg1], ['Prefit: Unfiltered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
			ax4.axes.xaxis.set_ticklabels([])

			ax5.errorbar(init_b3fr, init_b3re, yerr=init_b3Er, fmt='.', color='#D81B60', capsize=2)
			ax5.errorbar(init_b4fr, init_b4re, yerr=init_b4Er, fmt='.', color='green', capsize=2)
			ax5.grid()
			ax5.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b4_freq)+b3_bw/0.8))
			ax5.legend([leg1], ['Prefit: Filtered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
			ax5.axes.xaxis.set_ticklabels([])
			ax5.set_ylabel(r'ToA Residuals ($\mu$s)', fontweight='bold', fontsize=8)
		
			ax6.errorbar(final_b3fr, final_b3re, yerr=final_b3Er, fmt='.', color='#D81B60', capsize=2)
			ax6.errorbar(final_b4fr, final_b4re, yerr=final_b4Er, fmt='.', color='green', capsize=2)
			ax6.grid()
			ax6.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b4_freq)+b3_bw/0.8))
			ax6.legend([leg1], ['Postfit'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
			ax6.set_xlabel('Frequency (MHz)', fontweight='bold', fontsize=9)

			#ax7.errorbar(x1, y1, yerr=e1, fmt='.', color='#D81B60', capsize=2)
			#ax7.errorbar(x2, y2, yerr=e2, fmt='.', color='green', capsize=2)
			#ax7.errorbar(x3, y3, yerr=e3, fmt='.', color='#1E88E5', capsize=2)
			#ax7.grid()
			#ax7.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b4_freq)+b3_bw/0.8))
			#ax7.legend([leg1], ['Postfit'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
			#ax7.set_xlabel('Frequency (MHz)', fontweight='bold', fontsize=9)				
			#ax7.set_ylabel('P.fit.res($\mu$s))', fontweight='bold', fontsize=9)

			fig.suptitle('Source: PSR %s; MJD: %.4f; Prefit Wrms: %.2f $\mu$s; Postfit Wrms: %.2f $\mu$s\nMedian ToA Err: %.2f $\mu$s; DM: %.6f $\pm$ %.6f pc cm$^{-3}$;  Reduced $\chi^2$: %.2f'%(ar_psr, ar_mjd, prefit_rms[0], postfit_rms[0], med_toaE[0], dm[0], dmerror[0], chisq[0]), fontsize=10, fontweight='bold')


			dirplot=os.path.join(pwd,ar_psr+"_"+ar_tel+"_plots")
			if not os.path.exists(dirplot):
		   		os.makedirs(dirplot)
			plotfile=dirplot+"/"+ar_psr+"_"+str(ar_mjd)+"_"+str(centre_freqs[0])+"_"+ar_tel+"_DMfitResid.pdf"
			plt.savefig(plotfile, format='pdf')
			plt.close()



		if (2000. < cfs < 2300.): 

			prof2Db4 = []
			profb4 = []
			prof2Db5 = []
			profb5 = []
			b4_bw = []; b5_bw = []; b4_freq = []; b5_freq = []; b4_nbin = []; b5_nbin = []
			for i in range(len(finalarchives)):
				if (525. < finalarchives[i].get_centre_frequency() < 1000.):
					ar_nchan = finalarchives[i].get_nchan()
					b4_nbin  = finalarchives[i].get_nbin()
					b4_bw = finalarchives[i].get_bandwidth()
					b4_freq = finalarchives[i].get_centre_frequency()
					prof2Db4 = finalarchives[i].get_data()[:,0,:,:].flatten().reshape(ar_nchan,b4_nbin)
					prof = finalarchives[i].clone()
					prof.fscrunch()
					profb4 = prof.get_data().flatten()
					profb4 /= np.max(profb4)
				if (finalarchives[i].get_centre_frequency() > 1000.):
					ar_nchan = finalarchives[i].get_nchan()
					b5_nbin  = finalarchives[i].get_nbin()
					b5_bw = finalarchives[i].get_bandwidth()
					b5_freq = finalarchives[i].get_centre_frequency()
					prof2Db5 = finalarchives[i].get_data()[:,0,:,:].flatten().reshape(ar_nchan,b5_nbin)
					prof = finalarchives[i].clone()
					prof.fscrunch()
					profb5 = prof.get_data().flatten()
					profb5 /= np.max(profb5)
		
			condition = (orig_freqs < 1000.) & (orig_freqs > 525.)
			orig_b4fr = np.extract(condition,orig_freqs)
			orig_b4re = np.extract(condition,orig_resid)
			orig_b4Er = np.extract(condition,orig_toasE)

			condition = orig_freqs > 1000.
			orig_b5fr = np.extract(condition,orig_freqs)
			orig_b5re = np.extract(condition,orig_resid)
			orig_b5Er = np.extract(condition,orig_toasE)

			condition = (init_freqs < 1000.) & (init_freqs > 525.)
			init_b4fr = np.extract(condition,init_freqs)
			init_b4re = np.extract(condition,init_resid)
			init_b4Er = np.extract(condition,init_toasE)

			condition = init_freqs > 1000.
			init_b5fr = np.extract(condition,init_freqs)
			init_b5re = np.extract(condition,init_resid)
			init_b5Er = np.extract(condition,init_toasE)

			condition = (final_freqs < 1000.) & (final_freqs > 525.)
			final_b4fr = np.extract(condition,final_freqs)
			final_b4re = np.extract(condition,final_resid)
			final_b4Er = np.extract(condition,final_toasE)

			condition = final_freqs > 1000.
			final_b5fr = np.extract(condition,final_freqs)
			final_b5re = np.extract(condition,final_resid)
			final_b5Er = np.extract(condition,final_toasE)

			fig = plt.figure(3, figsize=(8, 6))
			fig.subplots_adjust(hspace=0.05)
			ax0 = plt.subplot2grid((9, 8), (0,0), rowspan=3, colspan=3)
			ax1 = plt.subplot2grid((9, 8), (3,0), rowspan=1, colspan=3)
			ax2 = plt.subplot2grid((9, 8), (5,0), rowspan=3, colspan=3)
			ax3 = plt.subplot2grid((9, 8), (8,0), rowspan=1, colspan=3)
		
			ax4 = plt.subplot2grid((9, 8), (0,4), colspan=4, rowspan=3)
			ax5 = plt.subplot2grid((9, 8), (3,4), colspan=4, rowspan=3)
			ax6 = plt.subplot2grid((9, 8), (6,4), colspan=4, rowspan=3)

			#ax7 = plt.subplot2grid((9, 8), (6,4), colspan=4, rowspan=3)

			leg1 = Rectangle((0, 0), 0, 0, alpha=0.0)
		
			ax0.imshow((np.sqrt(prof2Db5**2))**0.5, origin='lower', extent=(0,b5_nbin-1,(np.around(b5_freq)-b5_bw/2),(np.around(b5_freq)+b5_bw/2)), aspect='auto', cmap='hot')
			ax0.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=8)
			ax0.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
			ax1.plot(np.arange(b5_nbin, dtype=float),profb5, color='black', linewidth=0.5)
			ax1.set_xlim(0,b5_nbin-1)
			ax1.set_ylabel('Intensity', fontweight='bold', fontsize=8)

			ax2.imshow((np.sqrt(prof2Db4**2))**0.5, origin='lower', extent=(0,b4_nbin-1,(np.around(b4_freq)-b4_bw/2),(np.around(b4_freq)+b4_bw/2)), aspect='auto', cmap='hot')
			ax2.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=8)
			ax2.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
			ax3.plot(np.arange(b4_nbin, dtype=float),profb4, color='black', linewidth=0.5)
			ax3.set_xlim(0,b4_nbin-1)
			ax3.set_xlabel('Pulse Phase (bins)', fontweight='bold', fontsize=8)
			ax3.set_ylabel('Intensity', fontweight='bold', fontsize=8)

			ax4.errorbar(orig_b4fr, orig_b4re, yerr=orig_b4Er, fmt='.', color='green', capsize=2)
			ax4.errorbar(orig_b5fr, orig_b5re, yerr=orig_b5Er, fmt='.', color='#1E88E5', capsize=2)
			ax4.set_xlim((np.around(b4_freq)-b4_bw/0.8), (np.around(b5_freq)+b5_bw/0.8))
			ax4.grid()
			ax4.legend([leg1], ['Prefit: Unfiltered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
			ax4.axes.xaxis.set_ticklabels([])

			ax5.errorbar(init_b4fr, init_b4re, yerr=init_b4Er, fmt='.', color='green', capsize=2)
			ax5.errorbar(init_b5fr, init_b5re, yerr=init_b5Er, fmt='.', color='#1E88E5', capsize=2)
			ax5.grid()
			ax5.set_xlim((np.around(b4_freq)-b4_bw/0.8), (np.around(b5_freq)+b5_bw/0.8))
			ax5.legend([leg1], ['Prefit: Filtered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
			ax5.axes.xaxis.set_ticklabels([])
			ax5.set_ylabel(r'ToA Residuals ($\mu$s)', fontweight='bold', fontsize=8)
		
			ax6.errorbar(final_b4fr, final_b4re, yerr=final_b4Er, fmt='.', color='#D81B60', capsize=2)
			ax6.errorbar(final_b5fr, final_b5re, yerr=final_b5Er, fmt='.', color='#1E88E5', capsize=2)
			ax6.grid()
			ax6.set_xlim((np.around(b4_freq)-b4_bw/0.8), (np.around(b5_freq)+b5_bw/0.8))
			ax6.legend([leg1], ['Postfit'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
			ax6.set_xlabel('Frequency (MHz)', fontweight='bold', fontsize=9)

			#ax7.errorbar(x1, y1, yerr=e1, fmt='.', color='#D81B60', capsize=2)
			#ax7.errorbar(x2, y2, yerr=e2, fmt='.', color='green', capsize=2)
			#ax7.errorbar(x3, y3, yerr=e3, fmt='.', color='#1E88E5', capsize=2)
			#ax7.grid()
			#ax7.set_xlim((np.around(b4_freq)-b4_bw/0.8), (np.around(b5_freq)+b5_bw/0.8))
			#ax7.legend([leg1], ['Postfit'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
			#ax7.set_xlabel('Frequency (MHz)', fontweight='bold', fontsize=9)				
			#ax7.set_ylabel('P.fit.res($\mu$s))', fontweight='bold', fontsize=9)

			fig.suptitle('Source: PSR %s; MJD: %.4f; Prefit Wrms: %.2f $\mu$s; Postfit Wrms: %.2f $\mu$s\nMedian ToA Err: %.2f $\mu$s; DM: %.6f $\pm$ %.6f pc cm$^{-3}$;  Reduced $\chi^2$: %.2f'%(ar_psr, ar_mjd, prefit_rms[0], postfit_rms[0], med_toaE[0], dm[0], dmerror[0], chisq[0]), fontsize=10, fontweight='bold')


			dirplot=os.path.join(pwd,ar_psr+"_"+ar_tel+"_plots")
			if not os.path.exists(dirplot):
		   		os.makedirs(dirplot)
			plotfile=dirplot+"/"+ar_psr+"_"+str(ar_mjd)+"_"+str(centre_freqs[0])+"_"+ar_tel+"_DMfitResid.pdf"
			plt.savefig(plotfile, format='pdf')
			plt.close()

	if (len(finalarchives) > 2):
		prof2Db3 = []
		profb3 = []
		prof2Db5 = []
		profb5 = []
		prof2Db4 = []
		profb4 = []
		b3_bw = []; b5_bw = [];b4_bw = []; b3_freq = []; b5_freq = []; b4_freq = []; b3_nbin = []; b5_nbin = []; b4_nbin = []
		
		for i in range(len(finalarchives)):
			if (finalarchives[i].get_centre_frequency() < 500.):
				ar_nchan = finalarchives[i].get_nchan()
				b3_nbin  = finalarchives[i].get_nbin()
				b3_bw = finalarchives[i].get_bandwidth()
				b3_freq = finalarchives[i].get_centre_frequency()
				prof2Db3 = finalarchives[i].get_data()[:,0,:,:].flatten().reshape(ar_nchan,b3_nbin)
				prof = finalarchives[i].clone()
				prof.fscrunch()
				profb3 = prof.get_data().flatten()
				profb3 /= np.max(profb3)
			if (finalarchives[i].get_centre_frequency() > 1000.):
				ar_nchan = finalarchives[i].get_nchan()
				b5_nbin  = finalarchives[i].get_nbin()
				b5_bw = finalarchives[i].get_bandwidth()
				b5_freq = finalarchives[i].get_centre_frequency()
				prof2Db5 = finalarchives[i].get_data()[:,0,:,:].flatten().reshape(ar_nchan,b5_nbin)
				prof = finalarchives[i].clone()
				prof.fscrunch()
				profb5 = prof.get_data().flatten()
				profb5 /= np.max(profb5)
			if (525. < finalarchives[i].get_centre_frequency( )< 1000.):
				ar_nchan = finalarchives[i].get_nchan()
				b4_nbin  = finalarchives[i].get_nbin()
				b4_bw = finalarchives[i].get_bandwidth()
				b4_freq = finalarchives[i].get_centre_frequency()
				prof2Db4 = finalarchives[i].get_data()[:,0,:,:].flatten().reshape(ar_nchan,b4_nbin)
				prof = finalarchives[i].clone()
				prof.fscrunch()
				profb4 = prof.get_data().flatten()
				profb4 /= np.max(profb4)
		condition = orig_freqs < 500.
		orig_b3fr = np.extract(condition,orig_freqs)
		orig_b3re = np.extract(condition,orig_resid)
		orig_b3Er = np.extract(condition,orig_toasE)
        
		condition = (orig_freqs < 1000.) & (orig_freqs > 525.)
		orig_b4fr = np.extract(condition,orig_freqs)
		orig_b4re = np.extract(condition,orig_resid)
		orig_b4Er = np.extract(condition,orig_toasE)


		condition = orig_freqs > 1000.
		orig_b5fr = np.extract(condition,orig_freqs)
		orig_b5re = np.extract(condition,orig_resid)
		orig_b5Er = np.extract(condition,orig_toasE)

		condition = init_freqs < 500.
		init_b3fr = np.extract(condition,init_freqs)
		init_b3re = np.extract(condition,init_resid)
		init_b3Er = np.extract(condition,init_toasE)
        
		condition = (init_freqs < 1000.) & (init_freqs > 525.)
		init_b4fr = np.extract(condition,init_freqs)
		init_b4re = np.extract(condition,init_resid)
		init_b4Er = np.extract(condition,init_toasE)
        
		condition = init_freqs > 1000.
		init_b5fr = np.extract(condition,init_freqs)
		init_b5re = np.extract(condition,init_resid)
		init_b5Er = np.extract(condition,init_toasE)

		condition = final_freqs < 500.
		final_b3fr = np.extract(condition,final_freqs)
		final_b3re = np.extract(condition,final_resid)
		final_b3Er = np.extract(condition,final_toasE)

		condition = (final_freqs < 1000.) & (final_freqs > 525.)
		final_b4fr = np.extract(condition,final_freqs)
		final_b4re = np.extract(condition,final_resid)
		final_b4Er = np.extract(condition,final_toasE)

		condition = final_freqs > 1000.
		final_b5fr = np.extract(condition,final_freqs)
		final_b5re = np.extract(condition,final_resid)
		final_b5Er = np.extract(condition,final_toasE)

		fig = plt.figure(3, figsize=(8, 6))
		fig.subplots_adjust(hspace=0.05)
		ax0 = plt.subplot2grid((9, 8), (0,0), rowspan=1, colspan=3)
		ax1 = plt.subplot2grid((9, 8), (1,0), rowspan=1, colspan=3)
		axb40 = plt.subplot2grid((9, 8), (3,0), rowspan=1, colspan=3)
		axb41 = plt.subplot2grid((9, 8), (4,0), rowspan=1, colspan=3)
		ax2 = plt.subplot2grid((9, 8), (7,0), rowspan=1, colspan=3)
		ax3 = plt.subplot2grid((9, 8), (8,0), rowspan=1, colspan=3)
		
		ax4 = plt.subplot2grid((9, 8), (0,4), colspan=4, rowspan=3)
		ax5 = plt.subplot2grid((9, 8), (3,4), colspan=4, rowspan=3)
		ax6 = plt.subplot2grid((9, 8), (6,4), colspan=4, rowspan=3)

		#ax7 = plt.subplot2grid((9, 8), (6,4), colspan=4, rowspan=3)

		leg1 = Rectangle((0, 0), 0, 0, alpha=0.0)
		
		ax0.imshow((np.sqrt(prof2Db5**2))**0.5, origin='lower', extent=(0,b5_nbin-1,(np.around(b5_freq)-b5_bw/2),(np.around(b5_freq)+b5_bw/2)), aspect='auto', cmap='hot')
		ax0.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=6)
		ax0.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
		
		ax1.plot(np.arange(b5_nbin, dtype=float),profb5, color='black', linewidth=0.5)
		ax1.set_xlim(0,b5_nbin-1)
		ax1.set_ylabel('Intensity', fontweight='bold', fontsize=6)

		axb40.imshow((np.sqrt(prof2Db4**2))**0.5, origin='lower', extent=(0,b4_nbin-1,(np.around(b4_freq)-b4_bw/2),(np.around(b4_freq)+b4_bw/2)), aspect='auto', cmap='hot')
		axb40.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=6)
		axb40.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
		
		axb41.plot(np.arange(b4_nbin, dtype=float),profb4, color='black', linewidth=0.5)
		axb41.set_xlim(0,b4_nbin-1)
		axb41.set_ylabel('Intensity', fontweight='bold', fontsize=6)


		ax2.imshow((np.sqrt(prof2Db3**2))**0.5, origin='lower', extent=(0,b3_nbin-1,(np.around(b3_freq)-b3_bw/2),(np.around(b3_freq)+b3_bw/2)), aspect='auto', cmap='hot')
		ax2.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=6)
		ax2.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
		
		ax3.plot(np.arange(b3_nbin, dtype=float),profb3, color='black', linewidth=0.5)
		ax3.set_xlim(0,b3_nbin-1)
		ax3.set_xlabel('Pulse Phase (bins)', fontweight='bold', fontsize=6)
		ax3.set_ylabel('Intensity', fontweight='bold', fontsize=6)
				

		ax4.errorbar(orig_b3fr, orig_b3re, yerr=orig_b3Er, fmt='.', color='#D81B60', capsize=2)
		ax4.errorbar(orig_b4fr, orig_b4re, yerr=orig_b4Er, fmt='.', color='green', capsize=2)
		ax4.errorbar(orig_b5fr, orig_b5re, yerr=orig_b5Er, fmt='.', color='#1E88E5', capsize=2)
		ax4.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b5_freq)+b3_bw/0.8))
		ax4.grid()
		ax4.legend([leg1], ['Prefit: Unfiltered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		ax4.axes.xaxis.set_ticklabels([])

		ax5.errorbar(init_b3fr, init_b3re, yerr=init_b3Er, fmt='.', color='#D81B60', capsize=2)
		ax5.errorbar(init_b4fr, init_b4re, yerr=init_b4Er, fmt='.', color='green', capsize=2)
		ax5.errorbar(init_b5fr, init_b5re, yerr=init_b5Er, fmt='.', color='#1E88E5', capsize=2)
		ax5.grid()
		ax5.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b5_freq)+b3_bw/0.8))
		ax5.legend([leg1], ['Prefit: Filtered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		ax5.axes.xaxis.set_ticklabels([])
		ax5.set_ylabel(r'ToA Residuals ($\mu$s)', fontweight='bold', fontsize=8)
		
		ax6.errorbar(final_b3fr, final_b3re, yerr=final_b3Er, fmt='.', color='#D81B60', capsize=2)
		ax6.errorbar(final_b4fr, final_b4re, yerr=final_b4Er, fmt='.', color='green', capsize=2)
		ax6.errorbar(final_b5fr, final_b5re, yerr=final_b5Er, fmt='.', color='#1E88E5', capsize=2)
		ax6.grid()
		ax6.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b5_freq)+b3_bw/0.8))
		ax6.legend([leg1], ['Postfit'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		ax6.set_xlabel('Frequency (MHz)', fontweight='bold', fontsize=9)

		#ax7.errorbar(x1, y1, yerr=e1, fmt='.', color='#D81B60', capsize=2)
		#ax7.errorbar(x2, y2, yerr=e2, fmt='.', color='green', capsize=2)
		#ax7.errorbar(x3, y3, yerr=e3, fmt='.', color='#1E88E5', capsize=2)
		#ax7.grid()
		#ax7.set_xlim((np.around(b3_freq)-b3_bw/0.8), (np.around(b5_freq)+b3_bw/0.8))
		#ax7.legend([leg1], ['Postfit'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		#ax7.set_xlabel('Frequency (MHz)', fontweight='bold', fontsize=9)				
		#ax7.set_ylabel('P.fit.res($\mu$s))', fontweight='bold', fontsize=9)

		fig.suptitle('Source: PSR %s; MJD: %.4f; Prefit Wrms: %.2f $\mu$s; Postfit Wrms: %.2f $\mu$s\nMedian ToA Err: %.2f $\mu$s; DM: %.6f $\pm$ %.6f pc cm$^{-3}$;  Reduced $\chi^2$: %.2f'%(ar_psr, ar_mjd, prefit_rms[0], postfit_rms[0], med_toaE[0], dm[0], dmerror[0], chisq[0]), fontsize=10, fontweight='bold')


		dirplot=os.path.join(pwd,ar_psr+"_"+ar_tel+"_plots")
		if not os.path.exists(dirplot):
		   os.makedirs(dirplot)
		plotfile=dirplot+"/"+ar_psr+"_"+str(ar_mjd)+"_"+str(centre_freqs[0])+"_"+ar_tel+"_DMfitResid.pdf"
		plt.savefig(plotfile, format='pdf')
		plt.close()

	if (len(finalarchives) == 1):
		prof2D = []
		prof1D = []
		condition = orig_freqs < 500.
			
		ar_nchan = finalarchives[0].get_nchan()
		ar_nbin  = finalarchives[0].get_nbin()
		ar_bw = finalarchives[0].get_bandwidth()
		ar_freq = finalarchives[0].get_centre_frequency()
		prof2D = finalarchives[0].get_data()[:,0,:,:].flatten().reshape(ar_nchan,ar_nbin)
		prof = finalarchives[0].clone()
		prof.fscrunch()
		prof1D = prof.get_data().flatten()
		prof1D /= np.max(prof1D)

		fig,axs = plt.subplots(2, sharex=True, figsize=(8, 6))
		fig.subplots_adjust(hspace=0.05)
		ax0 = plt.subplot2grid((9, 8), (0,0), rowspan=7, colspan=3)
		ax1 = plt.subplot2grid((9, 8), (7,0), rowspan=2, colspan=3)

		ax2 = plt.subplot2grid((9, 8), (0,4), colspan=4, rowspan=3)
		ax3 = plt.subplot2grid((9, 8), (3,4), colspan=4, rowspan=3)
		ax4 = plt.subplot2grid((9, 8), (6,4), colspan=4, rowspan=3)
		leg1 = Rectangle((0, 0), 0, 0, alpha=0.0)
		ax0.imshow((np.sqrt(prof2D**2))**0.5, origin='lower', extent=(0,ar_nbin-1,(np.around(ar_freq)-ar_bw/2),(np.around(ar_freq)+ar_bw/2)), aspect='auto', cmap='hot')
		ax0.set_ylabel('Frequency (MHz)', fontweight='bold', fontsize=8)
		ax0.tick_params(axis='x', which='both', bottom=True, top=True, labelbottom=False)
		ax1.plot(np.arange(ar_nbin, dtype=float),prof1D, color='black', linewidth=0.5)
		ax1.set_xlim(0,ar_nbin-1)
		ax1.set_ylabel('Intensity', fontweight='bold', fontsize=8)
		ax1.set_xlabel('Pulse Phase (bins)', fontweight='bold', fontsize=8)


		if ar_freq < 500.:
			ax2.errorbar(orig_freqs, orig_resid, yerr=orig_toasE, fmt='.', color='#D81B60', capsize=2)
		if ar_freq > 1000.:
			ax2.errorbar(orig_freqs, orig_resid, yerr=orig_toasE, fmt='.', color='#1E88E5', capsize=2)
		if 525. < ar_freq < 1000.:
			ax2.errorbar(orig_freqs, orig_resid, yerr=orig_toasE, fmt='.', color='#1E88E5', capsize=2)
		ax2.set_xlim((np.around(ar_freq)-ar_bw/1.5), (np.around(ar_freq)+ar_bw/1.5))
		ax2.grid()
		ax2.legend([leg1], ['Prefit: Unfiltered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		ax2.axes.xaxis.set_ticklabels([])

		if ar_freq < 500.:
			ax3.errorbar(init_freqs, init_resid, yerr=init_toasE, fmt='.', color='#D81B60', capsize=2)
		if ar_freq > 1000.:
			ax3.errorbar(init_freqs, init_resid, yerr=init_toasE, fmt='.', color='#1E88E5', capsize=2)
		if 525. < ar_freq < 1000.:
			ax3.errorbar(init_freqs, init_resid, yerr=init_toasE, fmt='.', color='#1E88E5', capsize=2)
		ax3.grid()
		ax3.set_xlim((np.around(ar_freq)-ar_bw/1.5), (np.around(ar_freq)+ar_bw/1.5))
		ax3.legend([leg1], ['Prefit: Filtered'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		ax3.axes.xaxis.set_ticklabels([])
		ax3.set_ylabel(r'ToA Residuals ($\mu$s)', fontweight='bold', fontsize=8)
		
		if ar_freq < 500.:
			ax4.errorbar(final_freqs, final_resid, yerr=final_toasE, fmt='.', color='#D81B60', capsize=2)
		if ar_freq > 1000.:
			ax4.errorbar(final_freqs, final_resid, yerr=final_toasE, fmt='.', color='#1E88E5', capsize=2)
		if 525. < ar_freq < 1000.:
			ax4.errorbar(final_freqs, final_resid, yerr=final_toasE, fmt='.', color='green', capsize=2)    
        
		ax4.grid()
		ax4.set_xlim((np.around(ar_freq)-ar_bw/1.5), (np.around(ar_freq)+ar_bw/1.5))
		ax4.legend([leg1], ['Postfit'], handlelength=0, handletextpad=0, loc='upper right', fontsize=10)
		ax4.set_xlabel('Frequency (MHz)', fontweight='bold', fontsize=9)
		fig.suptitle('Source: PSR %s; MJD: %.4f; Prefit Wrms: %.2f $\mu$s; Postfit Wrms: %.2f $\mu$s\nMedian ToA Err: %.2f $\mu$s; DM: %.6f $\pm$ %.6f pc cm$^{-3}$;  Reduced $\chi^2$: %.2f'%(ar_psr, ar_mjd, prefit_rms[0], postfit_rms[0], med_toaE[0], dm[0], dmerror[0], chisq[0]), fontsize=10, fontweight='bold')
		
		#axs[0].set_ylabel(r'Prefit residuals ($\mu$s)', fontweight='bold', fontsize=12)
		#axs[0].errorbar(init_freqs,init_resid,init_toasE,fmt='.k', capsize=2)	
		#axs[1].errorbar(final_freqs,final_resid,final_toasE,fmt='.k', capsize=2)
		#axs[1].set_ylabel(r'Postfit residuals ($\mu$s)', fontweight='bold', fontsize=12)
		#axs[1].set_xlabel(r'Frequency (MHz)', fontweight='bold', fontsize=12)
		#fig.suptitle('Source: PSR %s; MJD: %.4f; Prefit Wrms: %.2f $\mu$s; Postfit Wrms: %.2f $\mu$s\nMedian ToA Err: %.2f $\mu$s; DM: %.6f $\pm$ %.6f pc cm$^{-3}$;  Reduced $\chi^2$: %.2f'%(ar_psr, ar_mjd, prefit_rms[0], postfit_rms[0], med_toaE[0], dm[0], dmerror[0], chisq[0]), fontsize=11, fontweight='bold')

		dirplot=os.path.join(pwd,ar_psr+"_"+ar_tel+"_plots")
		if not os.path.exists(dirplot):
		   os.makedirs(dirplot)
		plotfile=dirplot+"/"+ar_psr+"_"+str(ar_mjd)+"_"+str(centre_freqs[0])+"_"+ar_tel+"_DMfitResid.pdf"
		plt.savefig(plotfile, format='pdf')
		plt.close()
	
	if not quiet:
		print("Plots generated successfully.")
	import time
	end = time.time()
	total = end - start
	print ('\n-----------------------------------------------------------------------------')
	print ('MJD\t\tDM\t\tDMerr\t\tChisq\tC_Fr\tBW\tTel')
	print ('%.6f\t%.6f\t%.6f\t%.2f\t%.1f\t%.1f\t%s' % (ar_mjd, dm[0], dmerror[0], 
			chisq[0], centre_freqs[0], bw[0], ar_tel) )
	
	print ('-----------------------------------------------------------------------------')

	print("\nThe program took %.1f seconds to finish"%total)





# Function to get the final ToAs with proper uGMRT flags
def get_finalTOA(ar, std, ephemeris, templ_dm, quiet):
	
	if not quiet:
		print("Getting the final frequency resolved ToAs of %s..."%(ar.get_filename()))
	temppar = "tempnoefac.par"
	with open(ephemeris,"r") as ft:
		lines = ft.readlines()
	with open(temppar,"w") as ff:
		for line in lines:
			if not line.startswith("T2EFAC"):
				ff.write(line)
	ff.close()
	ar.set_ephemeris(temppar)
	ar.set_dispersion_measure(templ_dm)
	ar.update_model()
	ar_nchan = ar.get_nchan()
	snr = np.zeros(ar_nchan)
	artmp = ar.clone()
	for i in range(ar_nchan):
		snr[i] = artmp.get_Profile(0,0,i).snr()
	del (artmp)

	ar_psr = ar.get_source()
	ar_bw  = ar.get_bandwidth()
	ar_mjd = ar.get_Integration(0).get_start_time().in_days()
	ar_frq = ar.get_centre_frequency()
	pta_flag = '-pta InPTA'
	sys_flag = ''
	grp_flag = ''
	bnd_flag = ''
	ext_flag = ''
	cp = cdp_pa_mf(ar) 

	if (ar_mjd < 58600.):
		if (300. < ar_frq < 500.):
			grp_flag = '-group GM_GWB_'+cp[1]+'_'+str(ar_bw)+'_'+cp[0]+'_pre36'
			sys_flag = '-sys GM_GWB_'+cp[1]+'_'+str(ar_bw)+'_'+cp[0]
			bnd_flag = '-bandno 3'
		if (525. < ar_frq < 1000.):
			grp_flag = '-group GM_GWB_'+cp[1]+'_'+str(ar_bw)+'_'+cp[0]+'_pre36'
			sys_flag = '-sys GM_GWB_'+cp[1]+'_'+str(ar_bw)+'_'+cp[0]
			bnd_flag = '-bandno 4'
		if (1260. < ar_frq < 1460.):
			grp_flag = '-group GM_GWB_'+cp[1]+'_'+str(ar_bw)+'_'+cp[0]+'_pre36'
			sys_flag = '-sys GM_GWB_'+cp[1]+'_'+str(ar_bw)+'_'+cp[0]
			bnd_flag = '-bandno 5'
	if (ar_mjd > 58600.):
		if (300. < ar_frq < 500.):
			grp_flag = '-group GM_GWB_'+cp[1]+'_'+str(ar_bw)+'_'+cp[0]+'_post36'
			sys_flag = '-sys GM_GWB_'+cp[1]+'_'+str(ar_bw)+'_'+cp[0]
			bnd_flag = '-bandno 3'
		if (525. < ar_frq < 1000.):
			grp_flag = '-group GM_GWB_'+cp[1]+'_'+str(ar_bw)+'_'+cp[0]+'_post36'
			sys_flag = '-sys GM_GWB_'+cp[1]+'_'+str(ar_bw)+'_'+cp[0]
			bnd_flag = '-bandno 4'
		if (1260. < ar_frq < 1460.):
			grp_flag = '-group GM_GWB_'+cp[1]+'_'+str(ar_bw)+'_'+cp[0]+'_post36'
			sys_flag = '-sys GM_GWB_'+cp[1]+'_'+str(ar_bw)+'_'+cp[0]
			bnd_flag = '-bandno 5'

	if (ar_mjd > 58230.):
		ext_flag = '-cycle_post34'
	if (ar_mjd < 58230.):
		ext_flag = '-cycle_pre34'

	tmp_filename = os.path.basename(std.get_filename())
	std.set_filename(tmp_filename)
	ar_nbin  = ar.get_nbin()
	std_nbin = std.get_nbin()
	nbin = ar_nbin
	if (ar_nbin > std_nbin):
		nbin = std_nbin

	arrtim = psrchive.ArrivalTime()
	arrtim.set_shift_estimator('FDM:mcmc=1')
	arrtim.set_format('Tempo2')
	arrtim.set_format_flags('IPTA')
	std.bscrunch_to_nbin(nbin)
	std.pscrunch()
	std.tscrunch()
	std.dedisperse()
	arrtim.set_standard(std)
	tmp_filename = os.path.basename(ar.get_filename())
	ar.set_filename(tmp_filename)
	ar.bscrunch_to_nbin(nbin)
	ar.pscrunch()
	ar.tscrunch()
	arrtim.set_observation(ar)
	# Finding the ToAs and reading it into numpy arrays
	toas = arrtim.get_toas()
	toas_filtered = [x.split()[:] for x in toas]
	finaltoasfile = ar_psr+"_allToAs.tim"
	if not (os.path.isfile(finaltoasfile)):
		ft = open(finaltoasfile,"a")
		head="FORMAT 1"
		ft.write('%s\n' % head)
		ft.close()
	if (os.path.isfile(finaltoasfile)):
		ft = open(finaltoasfile,"a")

	for i in range(len(toas)):
		if (snr[i]>=8):
			ft.write('%s -prof_snr %.2f %s %s %s %s %s -chan %d\n' % (toas[i], snr[i], pta_flag, sys_flag, grp_flag, bnd_flag, ext_flag, i))
		if (snr[i]<8):
			ft.write('C %s -prof_snr %.2f %s %s %s %s %s -chan %d\n' % (toas[i], snr[i], pta_flag, sys_flag, grp_flag, bnd_flag, ext_flag, i))


	finaltoasfile1 = ar_psr+"_allToAs_copied.tim"
	if not (os.path.isfile(finaltoasfile1)):
		ft1 = open(finaltoasfile1,"a")
		head="FORMAT 1"
		ft1.write('%s\n' % head)
		ft1.close()
	if (os.path.isfile(finaltoasfile1)):
		ft1 = open(finaltoasfile1,"a")

	for i in range(len(toas)):
		if (snr[i]>=8):
			ft1.write('%s -prof_snr %.2f %s %s %s %s %s -chan %d\n' % (toas[i], snr[i], pta_flag, sys_flag, grp_flag, bnd_flag, ext_flag, i))
		if (snr[i]<8):
			ft1.write('C %s -prof_snr %.2f %s %s %s %s %s -chan %d\n' % (toas[i], snr[i], pta_flag, sys_flag, grp_flag, bnd_flag, ext_flag, i))

	if not quiet:
		print(" done!")




# Function to correct the backend delay for pre pinta V6.2 processed data
def Correct_delay(ar):
	if (ar.get_telescope() == 'GMRT' or ar.get_telescope() == 'gmrt'):
		ar_mjd = ar.get_Integration(0).get_start_time().in_days()
		ar_frq  = ar.get_centre_frequency()
		period = (ar.get_Integration(0).get_folding_period())
		# cycle 34-35
		if (ar_mjd >= 58230. and ar_mjd < 58550.):
			if (ar_frq> 300. and ar_frq < 500.):
				ar.unload('temp1.ar')
				tt=os.popen('psredit -c be:delay=-4.02653184 -m temp1.ar').read()
				ar = psrchive.Archive_load('temp1.ar')
				ar.update_model()
				os.remove('temp1.ar')
			if (ar_frq> 1200. and ar_frq < 1500.):
				ar.unload('temp2.ar')
				tt=os.popen('psredit -c be:delay=-4.02653184 -m temp2.ar').read()
				ar = psrchive.Archive_load('temp2.ar')
				ar.update_model()
				os.remove('temp2.ar')
		# cycle 37
		if (ar_mjd >= 58810. and ar_mjd < 58991.):
			if (ar_frq> 300. and ar_frq < 500.):
				ar.unload('temp1.ar')
				tt=os.popen('psredit -c be:delay=-2.01326592 -m temp1.ar').read()
				ar = psrchive.Archive_load('temp1.ar')
				ar.update_model()
				os.remove('temp1.ar')
			if (ar_frq> 1200. and ar_frq < 1500.):
				ar.unload('temp2.ar')
				tt=os.popen('psredit -c be:delay=-1.34217728 -m temp2.ar').read()
				ar = psrchive.Archive_load('temp2.ar')
				ar.update_model()
				os.remove('temp2.ar')
	return(ar)



def cdp_pa_mf(ar1):
	fn = ar1.get_filename()
	a = ''
	z = []
	command = f"psredit -c be:config "+fn
	result = subprocess.run(command, shell=True, stdout= subprocess.PIPE, text=True)
	if result.returncode == 0:
		output_lines = result.stdout.strip().split('\n')
		last_value = output_lines[-1].split('|')[-1].strip()
		first_value_line = output_lines[0]
		value = int(first_value_line.split('=')[1].split('|')[0])
	#print(last_value)
		if last_value == '0':
			a = 'b0'
		elif last_value == '1':
			a = 'b1'
		else:
			print("unexpected value:", last_value)
	else:
		print("Error running psredit command:")
		print(result.stderr)
	z = [a,str(value)]
	return(z)





def dm_estimate(psrname,ephemeris):
	timfile = psrname+"_allToAs_copied.tim"
	eph = ephemeris
	with open (timfile , 'r') as file:
		lines = file.readlines()
	modified_lines = []
	for l in lines:
		if l.startswith('C'):
			modified_lines.append(l[1:].lstrip())
		else:
			modified_lines.append(l)
	temp_tim = 'temptim.txt'
	with open (temp_tim, 'w') as file:
		file.writelines(modified_lines)


	filename, freq, toa, toaerr, telescope = [], [], [], [], []
	with open(temp_tim, 'r') as file:
		lin = file.readlines()
	
	for line in lin[1:]:
		columns = line.split()
		filename.append((columns[0]))
		freq.append(float(columns[1]))
		toa.append(float(columns[2]))
		toaerr.append(float(columns[3]))
		telescope.append((columns[4]))
	#return(columns)
	filename_arr = np.array(filename)
	freq_arr = np.array(freq)
	toa_arr = np.array(toa)
	toaerr_arr = np.array(toaerr)
	telescope_arr = np.array(telescope)

	#print(freq_arr)
	condition1 = toaerr_arr != 0.
	filenamenew = np.extract(condition1,filename_arr)
	freqnew = np.extract(condition1,freq_arr)
	toanew = np.extract(condition1,toa_arr)
	toaerrnew = np.extract(condition1,toaerr_arr)
	telescopenew = np.extract(condition1,telescope_arr)
	#print(freqnew)
##################################
	range_1_indicesN = np.where((freqnew >= 300) & (freqnew <= 505))
	range_2_indicesN = np.where((freqnew >= 550) & (freqnew <= 855))
	range_3_indicesN = np.where((freqnew >= 1160) & (freqnew <= 1465))
###	
	fr_b31 = freqnew[range_1_indicesN]
	filename_b31 = filenamenew[range_1_indicesN]
	toa_b31 = toanew[range_1_indicesN]
	toaerr_b31 = toaerrnew[range_1_indicesN]
	telescope_b31 = telescopenew[range_1_indicesN]

	conditionb3 = toaerr_b31 < 3*np.median(toaerr_b31)

	fr_b32 = np.extract(conditionb3,fr_b31)
	filename_b32 = np.extract(conditionb3,filename_b31)
	toa_b32 = np.extract(conditionb3,toa_b31)
	toaerr_b32 = np.extract(conditionb3,toaerr_b31)
	telescope_b32 = np.extract(conditionb3,telescope_b31)
	
###

###	
	fr_b41 = freqnew[range_2_indicesN]
	filename_b41 = filenamenew[range_2_indicesN]
	toa_b41 = toanew[range_2_indicesN]
	toaerr_b41 = toaerrnew[range_2_indicesN]
	telescope_b41 = telescopenew[range_2_indicesN]

	conditionb4 = toaerr_b41 < 3*np.median(toaerr_b41)

	fr_b42 = np.extract(conditionb4,fr_b41)
	filename_b42 = np.extract(conditionb4,filename_b41)
	toa_b42 = np.extract(conditionb4,toa_b41)
	toaerr_b42 = np.extract(conditionb4,toaerr_b41)
	telescope_b42 = np.extract(conditionb4,telescope_b41)
	
###

	fr_b51 = freqnew[range_3_indicesN]
	filename_b51 = filenamenew[range_3_indicesN]
	toa_b51 = toanew[range_3_indicesN]
	toaerr_b51 = toaerrnew[range_3_indicesN]
	telescope_b51 = telescopenew[range_3_indicesN]

	conditionb5 = toaerr_b51 < 3*np.median(toaerr_b51)

	fr_b52 = np.extract(conditionb5,fr_b51)
	filename_b52 = np.extract(conditionb5,filename_b51)
	toa_b52 = np.extract(conditionb5,toa_b51)
	toaerr_b52 = np.extract(conditionb5,toaerr_b51)
	telescope_b52 = np.extract(conditionb5,telescope_b51)

###########################################

	#condition2 = toaerrnew < 3*np.median(toaerrnew)
	filenamenew1 = np.concatenate((filename_b32, filename_b42, filename_b52))
	freqnew1 = np.concatenate((fr_b32, fr_b42, fr_b52))
	toanew1 = np.concatenate((toa_b32, toa_b42, toa_b52))
	toaerrnew1 = np.concatenate((toaerr_b32, toaerr_b42, toaerr_b52))
	telescopenew1 = np.concatenate((telescope_b32, telescope_b42, telescope_b52))
	#print(freqnew1)

	os.remove(temp_tim)

	tempfile1 = "tempnew.tim"

	f = open(tempfile1,"w+")
	head="FORMAT 1\n"
	f.write('%s' % head)
	for i in range(0,np.size(freqnew1)):
		f.write('%s %.8f %.18f %.6f %s \n' % (filenamenew1[i], freqnew1[i], 
				toanew1[i], toaerrnew1[i], telescopenew1[i]))

	f.close()

	





	tmp = os.popen("tempo2 -output general2 -f %s %s -s \"1111111 {freq} {pre} {err}\n\" | grep '1111111'" 
					% (ephemeris,tempfile1)).read()
	os.remove(tempfile1)

	# extracting the data from general2 output
	tmp1 = tmp.split('\n')
	

	#TErrtmp /= 1e+6
	freq = []
	residual =[]
	resid_err = []
	for entry in tmp1:
    		if entry:  # Check if the string is not empty
        		parts = entry.split()
        		freq.append(float(parts[1]))
			
	for entry in tmp1:
    		if entry:  # Check if the string is not empty
        		parts = entry.split()
        		residual.append(float(parts[2]))

	for entry in tmp1:
    		if entry:  # Check if the string is not empty
        		parts = entry.split()
        		resid_err.append(float(parts[3]))
	# Convert the list to a numpy array
	freq_array = np.array(freq)
	residual_array = np.array(residual)
	residerr_array = np.array(resid_err)

	#print("---for_check---")
	#print(freq_array,residual_array,residerr_array)
	#print("---for_check---")
	
	#residual_array = np.array(residual)
	residerr_array /= 1e+6
	#print(second_entries_array)

	range_1_indices = np.where((freq_array >= 300) & (freq_array <= 505))
	range_2_indices = np.where((freq_array >= 550) & (freq_array <= 855))
	range_3_indices = np.where((freq_array >= 1160) & (freq_array <= 1465))

	
	x1 = freq_array[range_1_indices]
	y1 = residual_array[range_1_indices]
	e1 = residerr_array[range_1_indices]
	toanew21 = toanew1[range_1_indices]
	toaerrnew21 = toaerrnew1[range_1_indices]
	filenamenew21 = filenamenew1[range_1_indices]
	telescopenew21 = telescopenew1[range_1_indices]

	x2 = freq_array[range_2_indices]
	y2 = residual_array[range_2_indices]
	e2 = residerr_array[range_2_indices]
	toanew22 = toanew1[range_2_indices]
	toaerrnew22 = toaerrnew1[range_2_indices]
	filenamenew22 = filenamenew1[range_2_indices]
	telescopenew22 = telescopenew1[range_2_indices]
	

	x3 = freq_array[range_3_indices]
	y3 = residual_array[range_3_indices]
	e3 = residerr_array[range_3_indices]
	toanew23 = toanew1[range_3_indices]
	toaerrnew23 = toaerrnew1[range_3_indices]
	filenamenew23 = filenamenew1[range_3_indices]
	telescopenew23 = telescopenew1[range_3_indices]

	xn1 = x1
	xn2 = x2
	xn3 = x3

	x1 = x1.reshape(-1,1)
	x2 = x2.reshape(-1,1)
	x3 = x3.reshape(-1,1)

	from sklearn import linear_model
	from sklearn.linear_model import HuberRegressor
	from sklearn.preprocessing import PolynomialFeatures
	from sklearn.pipeline import make_pipeline

	y1 *= 1e+6 
	y2 *= 1e+6
	y3 *= 1e+6

	#return(y1.size,y2.size,y3.size)
	print("Trying to remove outlier using Huber regression....")
	if (y1.size != 0 and  y2.size == 0 and y3.size == 0):
		toashift1 = (np.min(y1)*-1.5)


		y1 += toashift1
		

		e1 = e1*1e+6
		

		model1 = make_pipeline(PolynomialFeatures(2), HuberRegressor())
		

		try:
			model1.fit(x1,y1,huberregressor__sample_weight=np.ravel(1./e1))
			y_pred1 = model1.predict(x1)
			residuals1 = y1 - y_pred1
			median1 = np.median(residuals1)
			MAD1 = np.median(np.abs(residuals1-np.median(residuals1)))/0.6744897501960817
			condition1 = (residuals1 > median1 - 3*MAD1) & (residuals1 < median1 + 3*MAD1)
			
			freqf1 = np.extract(condition1,x1)		
			toaf1 = np.extract(condition1,toanew21)
			toaerrf1 = np.extract(condition1,toaerrnew21)	
			filenamef1 = np.extract(condition1,filenamenew21)
			telescopef1 = np.extract(condition1,telescopenew21)
			print("Successfully removed outliers using Huber regression in B3.")

		except ValueError as ve1:
			print(f"HuberRegressor failed with error: {ve1}")
			print("Moving without Huber Regression in B3.")
			freqf1 = xn1
			toaf1 = toanew21
			toaerrf1 = toaerrnew21
			filenamef1 = filenamenew21
			telescopef1 = telescopenew21

		
		filtered_array1 = [filenamef1, freqf1, toaf1, toaerrf1, telescopef1]
		filtered_array = np.array(filtered_array1)			
				
		
		dm3, dme3, chi3 = dmfit(filenamef1, freqf1, toaf1, toaerrf1, telescopef1, eph)
		
		dm = [dm3]
		dme = [dme3]
		chi = [chi3]
		
		#dminfo1 = [dm3, dme3, chi3]
		#dminfo = np.array(dminfo1)


		#return(xn1, xn3, freqf1, freqf3, toaf1, toaerrf1,filenamef1,telescopef1, np.size(freqf1),np.size(toaf1),dm3,dme3,chi3,dm5,dme5,chi5,dm3p5,dme3p5,chi3p5)

	if (y1.size == 0 and y2.size != 0 and y3.size == 0):
		toashift2 = (np.min(y2)*-1.5)


		y2 += toashift2
		

		e2 = e2*1e+6
		

		model2 = make_pipeline(PolynomialFeatures(2), HuberRegressor())
		

		try:
			model2.fit(x2,y2,huberregressor__sample_weight=np.ravel(1./e2))
			y_pred2 = model2.predict(x2)
			residuals2 = y2 - y_pred2
			median2 = np.median(residuals2)
			MAD2 = np.median(np.abs(residuals2-np.median(residuals2)))/0.6744897501960817
			condition2 = (residuals2 > median2 - 3*MAD2) & (residuals2 < median2 + 3*MAD2)

			freqf2 = np.extract(condition2,x2)
			toaf2 = np.extract(condition2,toanew22)
			toaerrf2 = np.extract(condition2,toaerrnew22)	
			filenamef2 = np.extract(condition2,filenamenew22)
			telescopef2 = np.extract(condition2,telescopenew22)
			print("Successfully removed outliers using Huber regression in B4.")
		except ValueError as ve2:
			print(f"HuberRegressor failed with error: {ve2}")
			print("Moving without Huber Regression in B4.")
			freqf2 = xn2
			toaf2 = toanew22
			toaerrf2 = toaerrnew22
			filenamef2 = filenamenew22
			telescopef2 = telescopenew22

		
			
		filtered_array1 = [filenamef2, freqf2, toaf2, toaerrf2, telescopef2]
		filtered_array = np.array(filtered_array1)		
		
		dm4, dme4, chi4 = dmfit(filenamef2, freqf2, toaf2, toaerrf2, telescopef2, eph)
		

		dm = [dm4]
		dme = [dme4]
		chi = [chi4]

		#dminfo1 = [dm4, dme4, chi4]
		#dminfo = np.array(dminfo1)


		#return(xn1, xn3, freqf1, freqf3, toaf1, toaerrf1,filenamef1,telescopef1, np.size(freqf1),np.size(toaf1),dm3,dme3,chi3,dm5,dme5,chi5,dm3p5,dme3p5,chi3p5)




	if (y1.size == 0 and y2.size == 0 and y3.size != 0):

		toashift3 = (np.min(y3)*-1.5)


		y3 += toashift3
		

		e3 = e3*1e+6
		

		model3 = make_pipeline(PolynomialFeatures(2), HuberRegressor())
		

		try:
			model3.fit(x3,y3,huberregressor__sample_weight=np.ravel(1./e3))
			y_pred3 = model3.predict(x3)
			residuals3 = y3 - y_pred3
			median3 = np.median(residuals3)
			MAD3 = np.median(np.abs(residuals3-np.median(residuals3)))/0.6744897501960817
			condition3 = (residuals3 > median3 - 3*MAD3) & (residuals3 < median3 + 3*MAD3)

			freqf3 = np.extract(condition3,x3)
			toaf3 = np.extract(condition3,toanew23)
			toaerrf3 = np.extract(condition3,toaerrnew23)	
			filenamef3 = np.extract(condition3,filenamenew23)
			telescopef3 = np.extract(condition3,telescopenew23)
			print("Successfully removed outliers using Huber regression in B5.")
		except ValueError as ve3:
			print(f"HuberRegressor failed with error: {ve3}")
			print("Moving without Huber Regression in B5.")
			freqf3 = xn3
			toaf3 = toanew23
			toaerrf3 = toaerrnew23
			filenamef3 = filenamenew23
			telescopef3 = telescopenew23
			
		filtered_array1 = [filenamef3, freqf3, toaf3, toaerrf3, telescopef3]
		filtered_array = np.array(filtered_array1)	
		
		dm5, dme5, chi5 = dmfit(filenamef3, freqf3, toaf3, toaerrf3, telescopef3, eph)
		


		dm = [dm5]
		dme = [dme5]
		chi = [chi5]

		#dminfo1 = [dm5, dme5, chi5]
		#dminfo = np.array(dminfo1)


		#return(xn1, xn3, freqf1, freqf3, toaf1, toaerrf1,filenamef1,telescopef1, np.size(freqf1),np.size(toaf1),dm3,dme3,chi3,dm5,dme5,chi5,dm3p5,dme3p5,chi3p5)






	if (y2.size == 0 and y1.size !=0 and y3.size !=0):
		
		toashift1 = (np.min(y1)*-1.5)
		toashift3 = (np.min(y3)*-1.5)

		y1 += toashift1
		y3 += toashift3

		e1 = e1*1e+6
		e3 = e3*1e+6

		model1 = make_pipeline(PolynomialFeatures(2), HuberRegressor())
		model3 = make_pipeline(PolynomialFeatures(2), HuberRegressor())

		try:
			model1.fit(x1,y1,huberregressor__sample_weight=np.ravel(1./e1))
			y_pred1 = model1.predict(x1)
			residuals1 = y1 - y_pred1
			median1 = np.median(residuals1)
			MAD1 = np.median(np.abs(residuals1-np.median(residuals1)))/0.6744897501960817
			condition1 = (residuals1 > median1 - 3*MAD1) & (residuals1 < median1 + 3*MAD1)
			
			freqf1 = np.extract(condition1,x1)		
			toaf1 = np.extract(condition1,toanew21)
			toaerrf1 = np.extract(condition1,toaerrnew21)	
			filenamef1 = np.extract(condition1,filenamenew21)
			telescopef1 = np.extract(condition1,telescopenew21)
			print("Successfully removed outliers using Huber regression in B3.")
		except ValueError as ve1:
			print(f"HuberRegressor failed with error: {ve1}")
			print("Moving without Huber Regression in B3.")
			freqf1 = xn1
			toaf1 = toanew21
			toaerrf1 = toaerrnew21
			filenamef1 = filenamenew21
			telescopef1 = telescopenew21

		try:
			model3.fit(x3,y3,huberregressor__sample_weight=np.ravel(1./e3))
			y_pred3 = model3.predict(x3)
			residuals3 = y3 - y_pred3
			median3 = np.median(residuals3)
			MAD3 = np.median(np.abs(residuals3-np.median(residuals3)))/0.6744897501960817
			condition3 = (residuals3 > median3 - 3*MAD3) & (residuals3 < median3 + 3*MAD3)

			freqf3 = np.extract(condition3,x3)
			toaf3 = np.extract(condition3,toanew23)
			toaerrf3 = np.extract(condition3,toaerrnew23)	
			filenamef3 = np.extract(condition3,filenamenew23)
			telescopef3 = np.extract(condition3,telescopenew23)
			print("Successfully removed outliers using Huber regression in B5.")
		except ValueError as ve3:
			print(f"HuberRegressor failed with error: {ve3}")
			print("Moving without Huber Regression in B5.")
			freqf3 = xn3
			toaf3 = toanew23
			toaerrf3 = toaerrnew23
			filenamef3 = filenamenew23
			telescopef3 = telescopenew23
		
		
		freq3p5 = np.concatenate((freqf1, freqf3))
		toaf3p5 = np.concatenate((toaf1,toaf3))
		toaerrf3p5 = np.concatenate((toaerrf1,toaerrf3))
		filenamef3p5 = np.concatenate((filenamef1,filenamef3))
		telescopef3p5 = np.concatenate((telescopef1,telescopef3))

		filtered_array1 = [filenamef3p5,freq3p5,toaf3p5,toaerrf3p5,telescopef3p5]
		filtered_array = np.array(filtered_array1)

		dm3, dme3, chi3 = dmfit(filenamef1, freqf1, toaf1, toaerrf1, telescopef1, eph)
		dm5, dme5, chi5 = dmfit(filenamef3, freqf3, toaf3, toaerrf3, telescopef3, eph)
		dm3p5, dme3p5, chi3p5 = dmfit(filenamef3p5,freq3p5,toaf3p5,toaerrf3p5,telescopef3p5,eph)

		
		dm = [dm3p5, dm3, dm5]
		dme = [dme3p5, dme3, dme5]
		chi = [chi3p5, chi3, chi5]
		#print(freqf1,freqf3,freq3p5)

		#prepost(filenamef3p5,freq3p5,toaf3p5,toaerrf3p5,telescopef3p5,eph,psrname)

		#dminfo1 = [dm3p5, dme3p5, chi3p5, dm3, dme3, chi3, dm5, dme5, chi5]
		#dminfo = np.array(dminfo1)


		#return(xn1, xn3, freqf1, freqf3, toaf1, toaerrf1,filenamef1,telescopef1, np.size(freqf1),np.size(toaf1),dm3,dme3,chi3,dm5,dme5,chi5,dm3p5,dme3p5,chi3p5)

	if (y3.size == 0 and y1.size !=0 and y2.size != 0):
		toashift1 = (np.min(y1)*-1.5)
		toashift2 = (np.min(y2)*-1.5)

		y1 += toashift1
		y2 += toashift2

		e1 = e1*1e+6
		e2 = e2*1e+6

		model1 = make_pipeline(PolynomialFeatures(2), HuberRegressor())
		model2 = make_pipeline(PolynomialFeatures(2), HuberRegressor())

		try:
			model1.fit(x1,y1,huberregressor__sample_weight=np.ravel(1./e1))
			y_pred1 = model1.predict(x1)
			residuals1 = y1 - y_pred1
			median1 = np.median(residuals1)
			MAD1 = np.median(np.abs(residuals1-np.median(residuals1)))/0.6744897501960817
			condition1 = (residuals1 > median1 - 3*MAD1) & (residuals1 < median1 + 3*MAD1)
			
			freqf1 = np.extract(condition1,x1)		
			toaf1 = np.extract(condition1,toanew21)
			toaerrf1 = np.extract(condition1,toaerrnew21)	
			filenamef1 = np.extract(condition1,filenamenew21)
			telescopef1 = np.extract(condition1,telescopenew21)
			print("Successfully removed outliers using Huber regression in B3.")
		except ValueError as ve1:
			print(f"HuberRegressor failed with error: {ve1}")
			print("Moving without Huber Regression in B3.")
			freqf1 = xn1
			toaf1 = toanew21
			toaerrf1 = toaerrnew21
			filenamef1 = filenamenew21
			telescopef1 = telescopenew21

		try:
			model2.fit(x2,y2,huberregressor__sample_weight=np.ravel(1./e2))
			y_pred2 = model2.predict(x2)
			residuals2 = y2 - y_pred2
			median2 = np.median(residuals2)
			MAD2 = np.median(np.abs(residuals2-np.median(residuals2)))/0.6744897501960817
			condition2 = (residuals2 > median2 - 3*MAD2) & (residuals2 < median2 + 3*MAD2)

			freqf2 = np.extract(condition2,x2)
			toaf2 = np.extract(condition2,toanew22)
			toaerrf2 = np.extract(condition2,toaerrnew22)	
			filenamef2 = np.extract(condition2,filenamenew22)
			telescopef2 = np.extract(condition2,telescopenew22)
			print("Successfully removed outliers using Huber regression in B4.")
		except ValueError as ve2:
			print(f"HuberRegressor failed with error: {ve2}")
			print("Moving without Huber Regression in B4.")
			freqf2 = xn2
			toaf2 = toanew22
			toaerrf2 = toaerrnew22
			filenamef2 = filenamenew22
			telescopef2 = telescopenew22

		
		freq3p4 = np.concatenate((freqf1, freqf2))
		toaf3p4 = np.concatenate((toaf1,toaf2))
		toaerrf3p4 = np.concatenate((toaerrf1,toaerrf2))
		filenamef3p4 = np.concatenate((filenamef1,filenamef2))
		telescopef3p4 = np.concatenate((telescopef1,telescopef2))

		filtered_array1 = [filenamef3p4,freq3p4,toaf3p4,toaerrf3p4,telescopef3p4]
		filtered_array = np.array(filtered_array1)

		dm3, dme3, chi3 = dmfit(filenamef1, freqf1, toaf1, toaerrf1, telescopef1, eph)
		dm4, dme4, chi4 = dmfit(filenamef2, freqf2, toaf2, toaerrf2, telescopef2, eph)
		dm3p4, dme3p4, chi3p4 = dmfit(filenamef3p4,freq3p4,toaf3p4,toaerrf3p4,telescopef3p4,eph)



		dm = [dm3p4, dm3, dm4]
		dme = [dme3p4, dme3, dme4]
		chi = [chi3p4, chi3, chi4]


	if (y1.size == 0 and y2.size !=0 and y3.size != 0):
		toashift2 = (np.min(y2)*-1.5)
		toashift3 = (np.min(y3)*-1.5)

		y2 += toashift2
		y3 += toashift3

		e3 = e3*1e+6
		e2 = e2*1e+6

		model3 = make_pipeline(PolynomialFeatures(2), HuberRegressor())
		model2 = make_pipeline(PolynomialFeatures(2), HuberRegressor())

		try:
			model3.fit(x3,y3,huberregressor__sample_weight=np.ravel(1./e3))
			y_pred3 = model3.predict(x3)
			residuals3 = y3 - y_pred3
			median3 = np.median(residuals3)
			MAD3 = np.median(np.abs(residuals3-np.median(residuals3)))/0.6744897501960817
			condition3 = (residuals3 > median3 - 3*MAD3) & (residuals3 < median3 + 3*MAD3)
			
			freqf3 = np.extract(condition3,x3)		
			toaf3 = np.extract(condition3,toanew23)
			toaerrf3 = np.extract(condition3,toaerrnew23)	
			filenamef3 = np.extract(condition3,filenamenew23)
			telescopef3 = np.extract(condition3,telescopenew23)
			print("Successfully removed outliers using Huber regression in B5.")
		except ValueError as ve3:
			print(f"HuberRegressor failed with error: {ve3}")
			print("Moving without Huber Regression in B5.")
			freqf3 = xn3
			toaf3 = toanew23
			toaerrf3 = toaerrnew23
			filenamef3 = filenamenew23
			telescopef3 = telescopenew23

		try:
			model2.fit(x2,y2,huberregressor__sample_weight=np.ravel(1./e2))
			y_pred2 = model2.predict(x2)
			residuals2 = y2 - y_pred2
			median2 = np.median(residuals2)
			MAD2 = np.median(np.abs(residuals2-np.median(residuals2)))/0.6744897501960817
			condition2 = (residuals2 > median2 - 3*MAD2) & (residuals2 < median2 + 3*MAD2)

			freqf2 = np.extract(condition2,x2)
			toaf2 = np.extract(condition2,toanew22)
			toaerrf2 = np.extract(condition2,toaerrnew22)	
			filenamef2 = np.extract(condition2,filenamenew22)
			telescopef2 = np.extract(condition2,telescopenew22)
			print("Successfully removed outliers using Huber regression in B4.")
		except ValueError as ve2:
			print(f"HuberRegressor failed with error: {ve2}")
			print("Moving without Huber Regression in B4.")
			freqf2 = xn2
			toaf2 = toanew22
			toaerrf2 = toaerrnew22
			filenamef2 = filenamenew22
			telescopef2 = telescopenew22

		
		freq4p5 = np.concatenate((freqf2, freqf3))
		toaf4p5 = np.concatenate((toaf2,toaf3))
		toaerrf4p5 = np.concatenate((toaerrf2,toaerrf3))
		filenamef4p5 = np.concatenate((filenamef2,filenamef3))
		telescopef4p5 = np.concatenate((telescopef2,telescopef3))

		filtered_array1 = [filenamef4p5,freq4p5,toaf4p5,toaerrf4p5,telescopef4p5]
		filtered_array = np.array(filtered_array1)

		dm5, dme5, chi5 = dmfit(filenamef3, freqf3, toaf3, toaerrf3, telescopef3, eph)
		dm4, dme4, chi4 = dmfit(filenamef2, freqf2, toaf2, toaerrf2, telescopef2, eph)
		dm4p5, dme4p5, chi4p5 = dmfit(filenamef4p5,freq4p5,toaf4p5,toaerrf4p5,telescopef4p5,eph)



		dm = [dm4p5, dm4, dm5]
		dme = [dme4p5, dme4, dme5]
		chi = [chi4p5, chi4, chi5]


		#prepost(filenamef3p4,freq3p4,toaf3p4,toaerrf3p4,telescopef3p4,eph,psrname)
		

		#dminfo1 = [dm3p4, dme3p4, chi3p4, dm3, dme3, chi3, dm4, dme4, chi4]
		#dminfo = np.array(dminfo1)


	
		#return(xn1, xn3, freqf1, freqf3, toaf1, toaerrf1,filenamef1,telescopef1, np.size(freqf1),np.size(toaf1),dm3,dme3,chi3,dm4,dme4,chi4,dm3p4,dme3p4,chi3p4)

	if (y1.size != 0 and y2.size != 0 and y3.size != 0):
		toashift1 = (np.min(y1)*-1.5)
		toashift2 = (np.min(y2)*-1.5)
		toashift3 = (np.min(y3)*-1.5)

		y1 += toashift1
		y2 += toashift2
		y3 += toashift3

		e1 = e1*1e+6
		e2 = e2*1e+6
		e3 = e3*1e+6

		model1 = make_pipeline(PolynomialFeatures(2), HuberRegressor())
		model2 = make_pipeline(PolynomialFeatures(2), HuberRegressor())
		model3 = make_pipeline(PolynomialFeatures(2), HuberRegressor())

		try:
			model1.fit(x1,y1,huberregressor__sample_weight=np.ravel(1./e1))
			y_pred1 = model1.predict(x1)
			residuals1 = y1 - y_pred1
			median1 = np.median(residuals1)
			MAD1 = np.median(np.abs(residuals1-np.median(residuals1)))/0.6744897501960817
			condition1 = (residuals1 > median1 - 3*MAD1) & (residuals1 < median1 + 3*MAD1)
			
			freqf1 = np.extract(condition1,x1)		
			toaf1 = np.extract(condition1,toanew21)
			toaerrf1 = np.extract(condition1,toaerrnew21)	
			filenamef1 = np.extract(condition1,filenamenew21)
			telescopef1 = np.extract(condition1,telescopenew21)
			print("Successfully removed outliers using Huber regression in B3.")
		except ValueError as ve1:
			print(f"HuberRegressor failed with error: {ve1}")
			print("Moving without Huber Regression in B3.")
			freqf1 = xn1
			toaf1 = toanew21
			toaerrf1 = toaerrnew21
			filenamef1 = filenamenew21
			telescopef1 = telescopenew21
			
		try:
			model2.fit(x2,y2,huberregressor__sample_weight=np.ravel(1./e2))
			y_pred2 = model2.predict(x2)
			residuals2 = y2 - y_pred2
			median2 = np.median(residuals2)
			MAD2 = np.median(np.abs(residuals2-np.median(residuals2)))/0.6744897501960817
			condition2 = (residuals2 > median2 - 3*MAD2) & (residuals2 < median2 + 3*MAD2)

			freqf2 = np.extract(condition2,x2)
			toaf2 = np.extract(condition2,toanew22)
			toaerrf2 = np.extract(condition2,toaerrnew22)	
			filenamef2 = np.extract(condition2,filenamenew22)
			telescopef2 = np.extract(condition2,telescopenew22)
			print("Successfully removed outliers using Huber regression in B4.")
		except ValueError as ve2:
			print(f"HuberRegressor failed with error: {ve2}")
			print("Moving without Huber Regression in B4.")
			freqf2 = xn2
			toaf2 = toanew22
			toaerrf2 = toaerrnew22
			filenamef2 = filenamenew22
			telescopef2 = telescopenew22

		try:
			model3.fit(x3,y3,huberregressor__sample_weight=np.ravel(1./e3))
			y_pred3 = model3.predict(x3)
			residuals3 = y3 - y_pred3
			median3 = np.median(residuals3)
			MAD3 = np.median(np.abs(residuals3-np.median(residuals3)))/0.6744897501960817
			condition3 = (residuals3 > median3 - 3*MAD3) & (residuals3 < median3 + 3*MAD3)

			freqf3 = np.extract(condition3,x3)
			toaf3 = np.extract(condition3,toanew23)
			toaerrf3 = np.extract(condition3,toaerrnew23)	
			filenamef3 = np.extract(condition3,filenamenew23)
			telescopef3 = np.extract(condition3,telescopenew23)
			print("Successfully removed outliers using Huber regression in B5.")
		except ValueError as ve3:
			print(f"HuberRegressor failed with error: {ve3}")
			print("Moving without Huber Regression in B5.")
			freqf3 = xn3
			toaf3 = toanew23
			toaerrf3 = toaerrnew23
			filenamef3 = filenamenew23
			telescopef3 = telescopenew23

		

		freqf3p4 = np.concatenate((freqf1,freqf2))
		toaf3p4 = np.concatenate((toaf1,toaf2))
		toaerrf3p4 = np.concatenate((toaerrf1,toaerrf2))
		filenamef3p4 = np.concatenate((filenamef1,filenamef2))
		telescopef3p4 = np.concatenate((telescopef1,telescopef2))

		freqf3p5 = np.concatenate((freqf1,freqf3))
		toaf3p5 = np.concatenate((toaf1,toaf3))
		toaerrf3p5 = np.concatenate((toaerrf1,toaerrf3))
		filenamef3p5 = np.concatenate((filenamef1,filenamef3))
		telescopef3p5 = np.concatenate((telescopef1,telescopef3))

		freqf4p5 = np.concatenate((freqf2,freqf3))
		toaf4p5 = np.concatenate((toaf2,toaf3))
		toaerrf4p5 = np.concatenate((toaerrf2,toaerrf3))
		filenamef4p5 = np.concatenate((filenamef2,filenamef3))
		telescopef4p5 = np.concatenate((telescopef2,telescopef3))


		freqf3p4p5 = np.concatenate((freqf1,freqf2, freqf3))
		toaf3p4p5 = np.concatenate((toaf1,toaf2, toaf3))
		toaerrf3p4p5 = np.concatenate((toaerrf1,toaerrf2, toaerrf3))
		filenamef3p4p5 = np.concatenate((filenamef1,filenamef2,filenamef3))
		telescopef3p4p5 = np.concatenate((telescopef1,telescopef2,telescopef3))

		filtered_array1 = [filenamef3p4p5,freqf3p4p5,toaf3p4p5,toaerrf3p4p5,telescopef3p4p5]
		filtered_array = np.array(filtered_array1)

		dm3, dme3, chi3 = dmfit(filenamef1, freqf1, toaf1, toaerrf1, telescopef1, eph)
		dm4, dme4, chi4 = dmfit(filenamef2, freqf2, toaf2, toaerrf2, telescopef2, eph)
		dm5, dme5, chi5 = dmfit(filenamef3, freqf3, toaf3, toaerrf3, telescopef3, eph)
		dm3p4, dme3p4, chi3p4 = dmfit(filenamef3p4, freqf3p4, toaf3p4, toaerrf3p4, telescopef3p4, eph)
		dm3p5, dme3p5, chi3p5 = dmfit(filenamef3p5, freqf3p5, toaf3p5, toaerrf3p5, telescopef3p5, eph)
		dm4p5, dme4p5, chi4p5 = dmfit(filenamef4p5, freqf4p5, toaf4p5, toaerrf4p5, telescopef4p5, eph)
		dm3p4p5, dme3p4p5, chi3p4p5 = dmfit(filenamef3p4p5,freqf3p4p5,toaf3p4p5,toaerrf3p4p5,telescopef3p4p5,eph)

		dm = [dm3p4p5,dm3p4,dm3p5,dm4p5,dm3,dm4,dm5]
		dme = [dme3p4p5,dme3p4,dme3p5,dme4p5,dme3,dme4,dme5]
		chi = [chi3p4p5,chi3p4,chi3p5,chi4p5,chi3,chi4,chi5]
	

		#dminfo1 = [dm3p4p5, dme3p4p5, chi3p4p5, dm3p4, dme3p4, chi3p4, dm3p5, dme3p5, chi3p5, dm4p5, dme4p5, chi4p5, dm3, dme3, chi3, dm4, dme4, chi4, dm5, dme5, chi5]
		#dminfo = np.array(dminfo1)

		#prepost(filenamef3p4p5,freqf3p4p5,toaf3p4p5,toaerrf3p4p5,telescopef3p4p5,eph,psrname)

	#return(x1,y1,e1,x2,y2,e2,x3,y3,e3)
		#return(dm3, dme3, chi3,dm4, dme4, chi4,dm5, dme5, chi5,dm3p4p5, dme3p4p5, chi3p4p5)

	return(dm,dme,chi,freq_array,residual_array,residerr_array,filtered_array)




def dmfit(filename, freqs, toas, toaE, tel, ephemeris):
	tempfile = filename[0]+"_"+str(freqs[0])+"_toas.txt"
	f = open(tempfile,"w+")
	head="FORMAT 1\n"
	f.write('%s' % head)
	for i in range(np.size(freqs)):
		f.write('%s %.8f %.18f %.6f %s\n' % (filename[i], freqs[i], toas[i], toaE[i], tel[i]))
	f.close()
	awkstr = "-nofit -fit dm | grep 'DM (cm^-3 pc)'| awk \'{print $5,$6}\'"
	dmstr=os.popen("tempo2 -f %s %s %s" % (ephemeris, tempfile, awkstr)).read()
	(dm, dmerr) = dmstr.split()
	dmval = float(dm)
	dmverr = float(dmerr)
	# doing the fit again to read the chisquare
	chisqstr=os.popen("tempo2 -f %s %s -nofit -fit dm | grep 'Fit Chisq'| awk \'{print $9}\'" % (ephemeris, tempfile)).read()
	fitchisq = float(chisqstr)
	#tempo2_output = os.popen("tempo2 -f %s %s -nofit -fit dm -output general2 -s \"11111 {freq} {pre} {post} {res} {err}e-6\n\" | grep '11111'" 
	#				% (ephemeris,tempfile)).read()
	#lines1 = tempo2_output.splitlines()
	#postfit_resid = filename[0]+"_resid.txt" 
	#with open(postfit_resid, "w") as f:
	#	f.write('\n'.join(lines1))
	os.remove(tempfile)
	return(dmval,dmverr,fitchisq)



main()
