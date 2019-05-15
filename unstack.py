#!/usr/bin/env python

import subprocess
import sys

filename = sys.argv[1]

subprocess.call(['e2proc2d.py',filename,'{0}.mrc'.format(filename.split('.')[0]),'--unstacking'])
