import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


file1 = sys.argv[1]
file2 = sys.argv[2]
file3 = sys.argv[3]

tmpfile = open(file3,"w+")
counter = 0
with open(file1) as fp:
    while True:
        line = fp.readline()
        if not line:
            break
        if line.strip().split()[0][0:3] != 'DMX' :
            tmpfile.write(line)
        else :
            counter = counter + 1
            if len(line.strip().split())==4:
                lineout = 'DMX_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '    ' +line.strip().split()[3] + '\n'
            elif len(line.strip().split())==3:
                lineout = 'DMX_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '    ' +line.strip().split()[2] + '\n'
            else :
                lineout = 'DMX_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '\n'
            tmpfile.write(lineout)
            line = fp.readline()
            if not line:
                break
            lineout = 'DMXEP_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '\n'
            tmpfile.write(lineout)
            line = fp.readline()
            if not line:
                break
            lineout = 'DMXR1_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '\n'
            tmpfile.write(lineout)
            line = fp.readline()
            if not line:
                break
            lineout = 'DMXR2_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '\n'
            tmpfile.write(lineout)
            line = fp.readline()
            if not line:
                break
            lineout = 'DMXF1_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '\n'
            tmpfile.write(lineout)
            line = fp.readline()
            if not line:
                break
            lineout = 'DMXF2_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '\n'
            tmpfile.write(lineout)

with open(file2) as fp:
    while True:
        line = fp.readline()
        if not line:
            break
        if line.strip().split()[0][0:3] == 'DMX' :
            counter = counter + 1
            if len(line.strip().split())==4:
                lineout = 'DMX_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '    ' +line.strip().split()[3] + '\n'
            elif len(line.strip().split())==3:
                lineout = 'DMX_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '    ' +line.strip().split()[2] + '\n'
            else :
                lineout = 'DMX_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '\n'
            tmpfile.write(lineout)
            line = fp.readline()
            if not line:
                break
            lineout = 'DMXEP_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '\n'
            tmpfile.write(lineout)
            line = fp.readline()
            if not line:
                break
            lineout = 'DMXR1_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '\n'
            tmpfile.write(lineout)
            line = fp.readline()
            if not line:
                break
            lineout = 'DMXR2_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '\n'
            tmpfile.write(lineout)
            line = fp.readline()
            if not line:
                break
            lineout = 'DMXF1_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '\n'
            tmpfile.write(lineout)
            line = fp.readline()
            if not line:
                break
            lineout = 'DMXF2_'+str(counter).zfill(4)+'     '+line.strip().split()[1] + '\n'
            tmpfile.write(lineout)

tmpfile.close()
