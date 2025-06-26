import numpy as np
import matplotlib.pyplot as plt
import psrchive
import time
import sys
import os


infile = sys.argv[1]

Deb = ['UD4', 'UD6']
Orig = ['Original']
Deb = ['UD4', 'UD6', 'UD8', 'UD10', 'UD12', 'UD14', 'UD16', 'UD18', 'UD20'] 
Spl = ['UB103', 'UB105', 'UB202', 'UB204', 'UB206', 'UB208', 'UB301', 'UB303', 'UB305', 'UB307', 'UB309']
Wavelt = Orig + Deb + Spl

wsnr=[]
nfile = 1;
cmd = f'psrstat -jDF -c "snr=pdmp" -c snr {infile} > temp.dat'
print('[CMD] :', cmd)
os.system(cmd)
f = open('temp.dat','r')
frow = f.readline()
f.close()
snr = frow.strip('\n\r').split()[1].split("=")[1]
print(snr)
wsnr.append(snr)
cmd = f'psrplot -Dwave.ps/ps -jDF -p flux -c below:r="SNR : {snr} \n Wvlt : None Original" {infile}'
print('[CMD] :', cmd)
os.system(cmd)
cmd = f'ps2pdf wave.ps'
print('[CMD] :', cmd)
os.system(cmd)
cmd = f'mv wave.pdf wavelet.pdf'
print('[CMD] :', cmd)
os.system(cmd)


for ielem in range(len(Deb)) :
    nfile = nfile + 1
    cmd = f'psrsmooth -W -t {Deb[ielem]} -e {Deb[ielem]} {infile}'
    print('[CMD] :', cmd)
    os.system(cmd)
    cmd = f'psrstat -jDF -c "snr=pdmp" -c snr {infile}.{Deb[ielem]} > temp.dat'
    print('[CMD] :', cmd)
    os.system(cmd)
    f = open('temp.dat','r')
    frow = f.readline()
    f.close()
    snr = frow.strip('\n\r').split()[1].split("=")[1]
    print(snr)
    wsnr.append(snr)
    cmd = f'psrplot -Dwave.ps/ps -jDF -p flux -c below:r="SNR : {snr} \n Wvlt : {Deb[ielem]}" {infile}.{Deb[ielem]}'
    print('[CMD] :', cmd)
    os.system(cmd)
    cmd = f'ps2pdf wave.ps'
    print('[CMD] :', cmd)
    os.system(cmd)
    cmd = f'pdfunite wave.pdf wavelet.pdf wavetemp.pdf'
    print('[CMD] :', cmd)
    os.system(cmd)
    cmd=f'mv wavetemp.pdf wavelet.pdf'
    print('[CMD] :', cmd)
    os.system(cmd)


for ielem in range(len(Spl)) :
    nfile = nfile + 1
    cmd = f'psrsmooth -W -t {Spl[ielem]} -e {Spl[ielem]} {infile}'
    print('[CMD] :', cmd)
    os.system(cmd)
    cmd = f'psrstat -jDF -c "snr=pdmp" -c snr {infile}.{Spl[ielem]} > temp.dat'
    print('[CMD] :', cmd)
    os.system(cmd)
    f = open('temp.dat','r')
    frow = f.readline()
    f.close()
    snr = frow.strip('\n\r').split()[1].split("=")[1]
    print(snr)
    wsnr.append(snr)
    cmd = f'psrplot -Dwave.ps/ps -jDF -p flux -c below:r="SNR : {snr} \n Wvlt : {Spl[ielem]}" {infile}.{Spl[ielem]}'
    print('[CMD] :', cmd)
    os.system(cmd)
    cmd = f'ps2pdf wave.ps'
    print('[CMD] :', cmd)
    os.system(cmd)
    cmd = f'pdfunite wave.pdf wavelet.pdf wavetemp.pdf'
    print('[CMD] :', cmd)
    os.system(cmd)
    cmd=f'mv wavetemp.pdf wavelet.pdf'
    print('[CMD] :', cmd)
    os.system(cmd)


f = open('wavelet.txt', 'w')
for i in range(nfile) :
    print(i)
    if i == 0 :
        cmd = f'Original  ||  SNR  : {wsnr[i]}\n'
    else :
        cmd = f'Wavelet : {Wavelt[i]}  ||  SNR  : {wsnr[i]}\n'
    f.write(cmd)
f.close()
    




    

