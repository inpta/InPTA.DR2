import numpy as np
import sys
import os

filename = sys.argv[1]
name = filename.split('/')[-1].split('.tim')[0]
outfile = f'{name}.tempflag.tim'

file = open(filename, 'r').readlines()
file = [f.split('\n')[0] for f in file]

newfile = open(outfile, 'w+')

tempflag1 = f'-tempflag1 1'
tempflag2 = f'-tempflag2 1'

for line in file:
    if 'FORMAT' not in line and 'MODE' not in line:
        linevals = line.split(' ')
        line = f''
        for entry in linevals[:-3]:
            line += f'{entry} '
            
        line += f'{linevals[-2]} {linevals[-1]} {linevals[-3]}'
        newline = f'{line} {tempflag1} {tempflag2}\n'
        newfile.write(newline)
    else:
        newfile.write(f'{line}\n')

newfile.close()


