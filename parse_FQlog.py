#!/usr/bin/env python

import sys

###---------function: read the star file get the header, labels, and data -------------#######
def read_starfile(f):
    inhead = True
    alldata = open(f,'r').readlines()
    labelsdic = {}
    data = []
    header = []
    count = 0
    labcount = 0
    for i in alldata:
        if '_rln' in i and '#' in i:
            labelsdic[i.split('#')[0]] = labcount
            labcount+=1
        if inhead == True:
            header.append(i.strip("\n"))
            if '_rln' in i and '#' in i and  '_rln' not in alldata[count+1] and '#' not in alldata[count+1]:
                inhead = False
        elif len(i.split())>=1:
            data.append(i.split())
        count +=1
    
    return(labelsdic,header,data)
#---------------------------------------------------------------------------------------------#

# get variables
try:
    (labels,header,data) = read_starfile(sys.argv[2])
    login = open(sys.argv[1]).readlines()
    threshold = float(sys.argv[3])
except:
    sys.exit("USAGE: parse_FQlog.py <FQ logfile> <original star file> <threshold>")

# get teh particles from FQ log above threshold
partlist = []
for i in login:
    line = i.split()
    if float(line[-1]) > threshold:
        part='{1:06d}@{0}.mrcs'.format(line[0].split('-')[0], int(line[0].split('-')[1].split('.')[0]))
        partlist.append(part)
        print(part)

# write star file
output = open('FQfilt.star','w')
output.write(header[0])
for i in header[1:]:
    output.write('\n{0}'.format(i))

for i in data:
    if i[labels['_rlnImageName ']] in partlist:
        output.write('\n{0}'.format('   '.join(i)))
