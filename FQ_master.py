#!/usr/bin/env python

import subprocess
import sys
import glob

# get the files
try:
    files = sys.argv[1:]
except:
    sys.exit('USAGE: FQ_master.py <files search string>\nfiles need to be in the proper path for relion and this script should be run in the relion working directory')

# hello world
print('{0} files'.format(len(files)))
for i in files:
    # unstack the file
    print('\nunstacking {0}'.format(i))
    workingdir = '/'.join(i.split('/')[0:-1])
    subprocess.call(['e2proc2d.py',i,'{0}.mrc'.format(i.split('.')[0]),'--unstacking'])
    print('working directory: {0}/'.format(workingdir))
    unstackedfiles = glob.glob('{0}/*.mrc'.format(workingdir))
    print('analysing particles')
    subprocess.call(['./fibril_quality.py','--i','{0}/'.format(workingdir),'--apix','1.065','--noimages'])
    print ('\ncleaning up {0} unstacked files'.format(len(unstackedfiles)))
    n=0
    for i in unstackedfiles:
        subprocess.call(['rm',i])
        n+=1
        if n%10 == 0:
            sys.stdout.write('.')
            sys.stdout.flush()