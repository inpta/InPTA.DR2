import numpy as np
import matplotlib.pyplot as plt
import psrchive
import time
import sys
import os

infile = sys.argv[1]
sbin = int(sys.argv[2])
ebin = int(sys.argv[3])
outfile = sys.argv[4]

centfreq = infile.split('_')[-1].split('.')[0]

cmd = f'pdv -R -A {infile} > temp.{centfreq}.dat'
print('[CMD] :', cmd)
os.system(cmd)

f = open(f'temp.{centfreq}.dat','r')
frow = f.readlines()
f.close()
frows=[l.strip('\n\r') for l in frow]
nlen = len(frows)

flds = frows[0].split()
nchan = int(flds[7])
nbin = int(flds[11])
psrname = flds[3]
mjd = int(float(frows[1].split()[1]))

if sbin >= ebin:
    print(f'The value of start_bin (= {sbin}) cannot be greater than or equal to the end_bin (= {ebin}) value!')
    print('The process will terminate!')
    sys.exit()    

elif sbin >= nbin:
    print(f'The value start_bin = {sbin} cannot be greater than or equal to the nbin (= {nbin}) value!')
    print('The process will terminate!')
    sys.exit()

elif ebin > nbin:
    print(f'The value end_bin = {ebin} cannot be greater than the nbin (= {nbin}) value!')
    print('The process will terminate!')
    sys.exit()

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
plt.title(f'{psrname} : {mjd} Bandshape', fontweight='bold')
plt.savefig(f'{psrname}.{mjd}.{centfreq}.bandshape.pdf')

weight = np.zeros(nchan)
for ichan in range(nchan) :
    if ichan !=0 and ichan != nchan-1:
        if bandshape[ichan] > 0.01*np.max(bandshape):
            weight[ichan] = np.max(bandshape)/bandshape[ichan]

plt.figure(figsize=(8, 5), dpi=150)      
plt.plot(xax, weight, lw=1, color='g')
plt.xlabel('Central Frequencies [MHz]', fontweight='bold')
plt.ylabel('Weight', fontweight='bold')
plt.title(f'{psrname} : {mjd} Channel Weights', fontweight='bold')
plt.savefig(f'{psrname}.{mjd}.{centfreq}.weights.pdf')

ar = psrchive.Archive_load(infile)
for ichan in range(nchan) :
    prof = ar.get_Integration(0).get_Profile(0, ichan)
    for ibin in range(nbin) :
        prof[ibin] = prof[ibin] * weight[ichan]
ar.unload(outfile)
