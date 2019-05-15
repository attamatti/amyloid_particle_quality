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
for i in files[0:1]:
    # unsatck the file
    print(i)
    workingdir = '/'.join(i.split('/')[0:-1])
    subprocess.call(['./unstack.py',i])
    print(workingdir)
    unstackedfiles = glob.glob('{0}/*.mrc'.format(workingdir))
    print ('{0} unstacked files'.format(len(unstackedfiles)))
    subprocess.call(['./fibril_quality.py','--i','{0}/'.format(workingdir),'--apix','1.065','--noimages'])