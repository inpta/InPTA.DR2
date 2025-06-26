import numpy as np
import math as mt
import statistics as stat
import psrchive
import sys
import os
import matplotlib
import pandas as pd
from statsmodels.stats import diagnostic
matplotlib.use('Agg')
import matplotlib.pyplot as plt 
import matplotlib.pylab as pylab
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
import time
from scipy.integrate import simpson

import functools
print = functools.partial(print, flush=True)

start = time.time()

input_fits = str(sys.argv[1])
nchanmin = int(sys.argv[2])
nchanmax = int(sys.argv[3])
dm_time_series = str(sys.argv[4])
start_bin = int(sys.argv[5])
end_bin = int(sys.argv[6])

def bandshape_equalise(infile, sbin, ebin, outfile, path):
    '''
        Function to equalise bandshape of profile across different channels
        
        InParams : Input FITS file, start and end phase bins in the off-pulse region,
                   name of the output FITS file (after equalisation), path to the storage directory
        Output : Bandshape-equalised FITS file
    '''
    
    print(f'\nPerforming bandshape equalisation considering ({sbin}, {ebin}) phase bin off-pulse region.')
    cmd = f'pdv -R -A {infile} > temp.dat'
    os.system(cmd)

    f = open('temp.dat','r')
    frow = f.readlines()
    f.close()
    frows=[l.strip('\n\r') for l in frow]
    nlen = len(frows)

    flds = frows[0].split()
    nchan = int(flds[7])
    nbin = int(flds[11])
    psrname = flds[3]
    mjd = frows[1].split()[1].split('.')[0]
    
    ielem = 1
    bandshape = np.zeros(nchan)
    xax = np.zeros(nchan)
    for ichan in range(nchan) :
        xax[ichan] = frows[ielem].split()[5]
        ielem = ielem + 1
        for ibin in range(nbin) :
            if ibin >= sbin and ibin <= ebin :
                bandshape[ichan] = bandshape[ichan] + float(frows[ielem].split()[3])
            ielem = ielem + 1
            if ielem >= nlen :
                break
        if ielem >= nlen :
            break

    plt.figure(figsize=(8, 5), dpi=150)      
    plt.plot(xax, bandshape, lw=1, color='r')
    plt.xlabel('Central Frequencies [MHz]', fontweight='bold')
    plt.ylabel('Amplitude', fontweight='bold')
    plt.title(f'{psrname} Bandshape for MJD {mjd}', fontweight='bold')
    plt.savefig(f'{path}/{psrname}.{mjd}.bandshape.pdf')

    weight = np.zeros(nchan)
    for ichan in range(nchan) :
        if ichan !=0 and ichan != nchan-1:
            if bandshape[ichan] > 0.01*np.max(bandshape):
                weight[ichan] = np.max(bandshape)/bandshape[ichan]

    plt.figure(figsize=(8, 5), dpi=150)      
    plt.plot(xax, weight, lw=1, color='g')
    plt.xlabel('Central Frequencies [MHz]', fontweight='bold')
    plt.ylabel('Weight', fontweight='bold')
    plt.title(f'{psrname} Channel Weights for MJD {mjd}', fontweight='bold')
    plt.savefig(f'{path}/{psrname}.{mjd}.weights.pdf')

    ar = psrchive.Archive_load(infile)
    for ichan in range(nchan) :
        prof = ar.get_Integration(0).get_Profile(0, ichan)
        for ibin in range(nbin) :
            prof[ibin] = prof[ibin] * weight[ichan]
            
    ar.unload(outfile)
    return outfile

def dm_dedisperse(fitsfile, dm_time_series):
    '''
        Function to dedisperse the PSRFITS file based on the DM of that epoch extracted from
        the DMcalc DM time series
        
        InParams : Input FITS file and DM time series file
        Output : Null
    '''
    data = np.genfromtxt(open(dm_time_series, 'r'), delimiter='')
    mjd_data = data[:, 0].astype(int)
    dm_data = data[:, 1]
    
    arch = psrchive.Archive_load(fitsfile)
    mjd_fitsfile = int(arch.get_Integration(0).get_start_time().in_days())

    index = np.where(mjd_data == mjd_fitsfile)[0][0]
    dm_fitsfile = dm_data[index]

    cmd = 'pam -m --update_dm '+str(dm_fitsfile)+' '+fitsfile
    os.system(cmd)

    cmd = 'pam -m -D '+fitsfile
    os.system(cmd)

if __name__ == '__main__' :
    
    threshold = np.round(np.arange(0.1, 1.01, 0.1), decimals=2)
    print('\n----------------------------------------------\n')
    
    fits = bandshape_equalise(input_fits, start_bin, end_bin, 'bandeq.fits', './')
    
    for thres in threshold :
        print('Setting threshold = %s'%thres)
        print('\n----------------------------------------------')
        
        thres_dir = 'threshold_'+str(thres)
        
        try : os.mkdir(thres_dir)
        except OSError : pass
    
        #Running while loop to do normalisation both with area and peak of the profile
        nchan = nchanmin
        while (nchan <= nchanmax):

            chan_dir = thres_dir+'/nchan_'+str(nchan)
            try : os.mkdir(chan_dir)
            except OSError : pass

            try : os.mkdir(chan_dir+'/Plots')
            except OSError : pass
            try : os.mkdir(chan_dir+'/Files')
            except OSError : pass
            try : os.mkdir(chan_dir+'/Files/Norm_area_method')
            except OSError : pass
            try : os.mkdir(chan_dir+'/Files/Norm_peak_method')
            except OSError : pass
            
            print(f"\nNormalising profile with area and peak, no. of channels = {nchan} \n")
            
            dm_dedisperse(fits, dm_time_series)
            
            arch = psrchive.Archive_load(fits)
            arch.tscrunch()
            arch.remove_baseline()
            arch.centre_max_bin()
            arch.fscrunch_to_nchan(nchan)
            
            plt.figure(figsize=(15,7))
            plt.subplot(1,2,1)
            ax1 = plt.gca()
            ax1.set_title("Normalization by peak")
            plt.subplot(1,2,2)
            ax2 = plt.gca()
            ax2.set_title("Normalization by area")
            
            profdiffpeak=np.array([0]*(nchan-1),dtype=object)
            profdiffarea=np.array([0]*(nchan-1),dtype=object)
    
            for i in range(0,nchan-1):
                
                prof1 = arch.get_Profile(0,0,i).get_amps()
                maxprof1 = np.max(prof1)
                profpeaknorm1 = prof1/(maxprof1)  #normalization by peak
                
                # sumprof1 = sum(prof1)
                nbin1 = arch.get_Profile(0,0,i).get_nbin()
                # areanorm1 = sumprof1*nbin1
                # profareanorm1 = prof1/(areanorm1) #normalization by area

                # Using integral function
                # -----
                phase1 = np.arange(0, nbin1, 1)
                intg1 = simpson(prof1, x=phase1, dx=1)
                profareanorm1 = prof1/intg1
                # -----

                prof2 = arch.get_Profile(0,0,i+1).get_amps()
                maxprof2 = np.max(prof2)                
                profpeaknorm2 = prof2/(maxprof2)  #normalization by peak
                
                # sumprof2 = sum(prof2)
                nbin2 = arch.get_Profile(0,0,i+1).get_nbin()
                # areanorm2 = sumprof2*nbin2
                # profareanorm2 = prof2/(areanorm2) #normalization by area

                # Using integral function
                # -----
                phase2 = np.arange(0, nbin2, 1)
                intg2 = simpson(prof2, x=phase2, dx=1)
                profareanorm2 = prof2/intg2
                # -----
                
                profdiffarea[i] = profareanorm1 - profareanorm2
                profdiffpeak[i] = profpeaknorm1 - profpeaknorm2
                
                df_diff_area = pd.DataFrame(profdiffarea[i])
                df_diff_peak = pd.DataFrame(profdiffpeak[i])
                
                df_diff_area.to_csv(chan_dir+f"/Files/Norm_area_method/Profile_subbnd_diff_{i}_{i+1}.txt", index=False, sep=' ')
                df_diff_peak.to_csv(chan_dir+f"/Files/Norm_peak_method/Profile_subbnd_diff_{i}_{i+1}.txt", index=False, sep=' ')
    
            rmsdiffpeak=np.array([0]*(nchan-1),dtype=object)
            rmsdiffarea=np.array([0]*(nchan-1),dtype=object)
    
            for i in range(0,nchan-1):
                rmsdiffpeak[i]=mt.sqrt(stat.mean(profdiffpeak[i]**2))
                rmsdiffarea[i]=mt.sqrt(stat.mean(profdiffarea[i]**2))
                
            medrmsdiffpeak=stat.median(rmsdiffpeak)
            medrmsdiffarea=stat.median(rmsdiffarea)
            meanrmsdiffpeak=stat.mean(rmsdiffpeak)
            meanrmsdiffarea=stat.mean(rmsdiffarea)
            
            print(f"Median RMS residual (peak norm) {medrmsdiffpeak}; Median RMS residual (area norm) {medrmsdiffarea} \n")
            
            # Although called madrms the quantity below is rms estimated from MAD
            madrmsdiffpeak=stat.median(abs(rmsdiffpeak-medrmsdiffpeak))*1.4826
            madrmsdiffarea=stat.median(abs(rmsdiffarea-medrmsdiffarea))*1.4826
            
            print(f"MAD RMS residual (peak norm) {madrmsdiffpeak}; MAD RMS residual (area norm) {madrmsdiffarea} \n")
    
            for i in range(0, nchan-1):
                if (rmsdiffpeak[i]-meanrmsdiffpeak <= thres*madrmsdiffpeak):
                    ax1.plot(profdiffpeak[i], label = f"Subband {i} - {i+1}")
                    
                if (rmsdiffarea[i]-meanrmsdiffarea <= thres*madrmsdiffarea):
                    ax2.plot(profdiffarea[i], label = f"Subband {i} - {i+1}")
                    
                if ((rmsdiffpeak[i]-meanrmsdiffpeak <= thres*madrmsdiffpeak) and (rmsdiffarea[i]-meanrmsdiffarea <= thres*madrmsdiffarea)):
                    aarea, pvaluearea = diagnostic.normal_ad(profdiffarea[i])
                    print(f"\tP-value: {pvaluearea}, area norm, sub {i}-{i+1}")
                    apeak, pvaluepeak = diagnostic.normal_ad(profdiffpeak[i])
                    print(f"\tP-value: {pvaluepeak}, peak norm, sub {i}-{i+1}")
                    
                if ((rmsdiffpeak[i]-meanrmsdiffpeak <= thres*madrmsdiffpeak) and (rmsdiffarea[i]-meanrmsdiffarea > thres*madrmsdiffarea)):
                    apeak, pvaluepeak = diagnostic.normal_ad(profdiffpeak[i])
                    print(f"\tP-value: {pvaluepeak}, peak norm, sub {i}-{i+1}")
                    
                if ((rmsdiffpeak[i]-meanrmsdiffpeak > thres*madrmsdiffpeak) and (rmsdiffarea[i]-meanrmsdiffarea <= thres*madrmsdiffarea)):
                    aarea, pvaluearea = diagnostic.normal_ad(profdiffarea[i])
                    print(f"\tP-value: {pvaluearea}, area norm, sub {i}-{i+1}")
    
            ax2.legend(fontsize=8, bbox_to_anchor=(1.255555, 1), loc='upper right')
            plt.savefig(chan_dir+f"/Plots/Prof_Diff_nchan{nchan}.png")
            
            plt.figure(figsize=(18,9))
            plt.subplot(1,2,1)
            ax3 = plt.gca()
            ax3.set_title("Normalization by peak")
            plt.subplot(1,2,2)
            ax4 = plt.gca()
            ax4.set_title("Normalization by area")
    
            for i in range(0, nchan):
                prof = arch.get_Profile(0,0,i).get_amps()
                nbin = arch.get_Profile(0,0,i).get_nbin()
                
                # normarea = sum(prof)*nbin
                # profnormarea = prof/normarea

                # Using integral function
                # -----
                phase = np.arange(0, nbin, 1)
                intg = simpson(prof, x=phase, dx=1)
                profnormarea = prof/intg
                # -----
                
                normpeak = np.max(prof)
                profnormpeak = prof/normpeak
                
                np.savetxt(chan_dir+f"/Files/Norm_area_method/Norm_profile_subbnd{i}.txt", profnormarea, fmt='%s')
                np.savetxt(chan_dir+f"/Files/Norm_peak_method/Norm_profile_subbnd{i}.txt", profnormpeak, fmt='%s')
                
                ax3.plot(profnormpeak, label = f"Subband {i}")
                ax4.plot(profnormarea, label = f"Subband {i}")
    
            ax4.legend(fontsize=8, bbox_to_anchor=(1.2, 1), loc='upper right')
            plt.savefig(chan_dir+f"/Plots/Norm_Profile_nchan{nchan}.png")

            # Plotting the 2D normalised profiles (frequency-phase plots)
            # -------
            nbin = arch.get_Profile(0, 0, 0).get_nbin()
            chans = np.arange(0, nchan, 1)
            cen_freqs = np.round(np.linspace(500, 300, nchan), decimals=1)
            data_area = np.zeros((nbin, nchan))
            data_peak = np.zeros((nbin, nchan))

            for i in range(nchan):
                data_area[:, i] = np.genfromtxt(open(chan_dir+f'/Files/Norm_area_method/Norm_profile_subbnd{i}.txt', 'r'))
                data_peak[:, i] = np.genfromtxt(open(chan_dir+f'/Files/Norm_peak_method/Norm_profile_subbnd{i}.txt', 'r'))

            plt.figure(figsize=(10, 8), dpi=100)
            plt.imshow(data_area[:, ::-1].T, aspect='auto')
            plt.yticks(chans, labels=cen_freqs)
            plt.xlabel('Phase Bins', fontweight='bold', fontsize=13)
            plt.ylabel('Central Frequency (MHz)', fontweight='bold', fontsize=13)
            plt.title(f'2D profile (Area Normalised) for nchan = {nchan}', fontweight='bold', fontsize=14)
            plt.colorbar(label='Signal Strength')
            plt.savefig(chan_dir+f'/Plots/Area_norm_2Dprofile.pdf')

            plt.figure(figsize=(10, 8), dpi=100)
            plt.imshow(data_peak[:, ::-1].T, aspect='auto')
            plt.yticks(chans, labels=cen_freqs)
            plt.xlabel('Phase Bins', fontweight='bold', fontsize=13)
            plt.ylabel('Central Frequency (MHz)', fontweight='bold', fontsize=13)
            plt.title(f'2D profile (Peak Normalised) for nchan = {nchan}', fontweight='bold', fontsize=14)
            plt.colorbar(label='Signal Strength')
            plt.savefig(chan_dir+f'/Plots/Peak_norm_2Dprofile.pdf')
            # -------
            
            nchan = 2*nchan
    
        print('\n----------------------------------------------\n')

end = time.time()
print('The process took %ds to finish'%(end-start))
print('\n----------------------------------------------\n')







