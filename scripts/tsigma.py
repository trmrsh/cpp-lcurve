#!/usr/bin/env python

usage = \
"""
This script calculates a lower limit to the uncertainty in t0 that would be 
measured using a given model applied to a given light curve. Although a lower
limit, I have found this gives a pretty good estimate of the actual values that
is obtained in practice.

It works by calculating the derivative with respect to time of the light curve.
There is then an equation using these and the data uncertainties which gives an
estimate of the timing accuracy. 

Once it has finished it makes a plot of the model light curve and the derivatives 
(scaled so have similar order of magnitudes although the scaling can be adjusted
using the -d option). This can be useful in assessing the degree of numerical noise
your derivatives might be suffering from.
"""

import os
import argparse
import subprocess
import pylab as plt
import numpy as np

# Light curve command: adjust this to point to wherever you have installed
# lroche (try 'which lroche')

lroche = os.path.join(os.environ['TRM_SOFTWARE'],  'bin', 'lcurve', 'lroche')

parser = argparse.ArgumentParser(description=usage)

# positional
parser.add_argument('model', type=argparse.FileType('r'), help='light curve model file')
parser.add_argument('data',  type=argparse.FileType('r'), help='light curve data file')
parser.add_argument('delta', type=float, help='perturbation step in time used to calculate numerical derivatives')

parser.add_argument('-d', dest='dmax', type=float, default=1, 
                    help='Normalisation factor to use when scaling derivatives for plot')

# OK, done with arguments.
args = parser.parse_args()

model  = args.model.name
data   = args.data.name

args.model.close()
args.data.close()

# compute lightcurve to set the scaling factor
print('\nComputing model light curve ...\n')

subprocess.call([lroche,model,data,'0','12345','1','junk.dat','/null','yes'])

# read result
mod = np.loadtxt('junk.dat')

# read in the model in order to modify it
with open(model) as fp:
    lines = fp.readlines()

# look for the t0 line
for i, line in enumerate(lines):
    if line.startswith('t0'):
        tindex = i
        break
else:
    print('Failed to find t0 line')
    exit(1)

arr = lines[tindex].split()

# perturb t0 by adding delta to it
print('\nPerturbing t0 to t0+delta ...\n')
lines[tindex] = '%s = %17.11f %s %s %s\n' % (arr[0],float(arr[2])+args.delta,arr[3],arr[4],arr[5])

# write out
with open('junk.mod','w') as fp:
    fp.writelines(lines)

# compute noiseless light curve with same scaling factor as just set
subprocess.call([lroche,'junk.mod',data,'0','12345','1','junk.dat','/null','no','\\'])

# read in result
lcp = np.loadtxt('junk.dat')

# perturb by subtracting delta
print('\nPerturbing t0 to t0-delta ...\n')

lines[tindex] = '%s = %17.11f %s %s %s\n' % (arr[0],float(arr[2])-args.delta,arr[3],arr[4],arr[5])

with open('junk.mod','w') as fp:
    fp.writelines(lines)

subprocess.call([lroche,'junk.mod',data,'0','12345','1','junk.dat','/null','no','\\'])

lcm = np.loadtxt('junk.dat')

# compute derivative
deriv = (lcm[:,2]-lcp[:,2])/(2.*args.delta)

# load data
lc = np.loadtxt(data)

# rcompute chi**2
rchisq = (((lc[:,2]-mod[:,2])/lc[:,3])**2).sum()/(len(lc[:,2])-1)
print('Reduced chisq =',rchisq)

# compute lower limit to uncertainty on t0
sigma = 1./np.sqrt(((deriv/lc[:,3])**2).sum())
print('\n>> Estimated minimum uncertainty in t0  =',sigma)
print('>> and after scaling reduced chisq to 1 =',sigma*np.sqrt(rchisq))


deriv /= np.max(np.abs(deriv))*args.dmax
flux   = mod[:,2]/mod[:,2].mean()

# make a plot of the derivative and light curve (normalised to have a max abs value of 1)
t     = lcm[:,0]
plt.plot(t,flux)
plt.plot([t[0],t[-1]],[0,0],'k--')
plt.plot(t,deriv)
plt.plot(t,deriv**2)
plt.show()

