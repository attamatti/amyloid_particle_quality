#!/usr/bin/env python
#
# find micrographs with enhanced signal at 4.7-4.9 A resolution
# calculate radialy averaged power spectrum of image
# caclulate fit based on 6.5 - 5.0 A and 3.0 - 2.8 A   



import glob
import sys
import os
import matplotlib.pyplot as plt
from scipy import fftpack
from scipy import misc
import numpy as np
import struct
from numpy import *

# scipy is triggering a future warning
import warnings
warnings.filterwarnings("ignore")

vers = '0.1'

# for troubleshooting arrays
#np.set_printoptions(threshold=np.nan)


#------------- READ MRCs--------------
# from David Stokes - modified by MGI
class mrc_image:
    def __init__(self,filename):
        self.numbytes1=56           # 56 long ints
        self.numbytes2=80*10          # 10 header lines of 80 chars each
        self.filename=filename

def read(self):
    global imgs_dic
    input_image=open(self.filename,'rb')
    self.header1=input_image.read(self.numbytes1*4)
    self.header2=input_image.read(self.numbytes2)
    byte_pattern='=' + 'l'*self.numbytes1   #'=' required to get machine independent standard size
    self.dim=struct.unpack(byte_pattern,self.header1)[:3]   #(dimx,dimy,dimz)
    self.imagetype=struct.unpack(byte_pattern,self.header1)[3]  #0: 8-bit signed, 1:16-bit signed, 2: 32-bit float, 6: unsigned 16-bit (non-std)
    if (self.imagetype == 0):
        imtype='b'
    elif (self.imagetype == 1):
        imtype='h'
    elif (self.imagetype == 2):
        imtype='f4'
    elif (self.imagetype == 6):
        imtype='H'
    else:
        imtype='unknown'   #should put a fail here
        print('WARNING: MRC type unrecognized - check results carefully')
    if stack == False:
        input_image_dimension=(self.dim[(1)],self.dim[(0)])  #2D images assumed
        self.image_data=fromfile(file=input_image,dtype=imtype,count=self.dim[(0)]*self.dim[(1)]).reshape(input_image_dimension)
    elif stack == True:
        input_image_dimension=(self.dim[(2)],self.dim[(1)],self.dim[(0)])  #3D images assumed
        self.image_data=fromfile(file=input_image,dtype=imtype,count=self.dim[(0)]*self.dim[(1)]*self.dim[(2)]).reshape(input_image_dimension)
        #print('dimensions: {0} x {1} x {2}'.format(self.dim[(1)],self.dim[(0)],self.dim[(2)]))
        imgs_dic[self.filename] = self.dim[(2)]
        return(np.mean(self.image_data,axis=0))
    input_image.close()
    #print('dimensions: {0} x {1} x {2}'.format(self.dim[(1)],self.dim[(0)],self.dim[(2)]))
    return(self.image_data)
#---------------------------

#-----------------------------------------------------------------------------------#
errormsg = """USAGE: fibril_quality.py --i <path to directory containing micrographs/particle stacks> --apix <angstroms per pixel>
If the input images are stacks add '--stack' to the input line

Make sure to include the trailing / on the directory name; the program adds *.mrc to the end of it to find both .mrc or .mrcs files
The script operates on .mrc single images (micrographs) or .mrcs stacks (relion particle stacks) but not both at the same time
"""

class Arg(object):
    _registry = []
    def __init__(self, flag, value, req):
        self._registry.append(self)
        self.flag = flag
        self.value = value
        self.req = req

def make_arg(flag, value, req):
    Argument = Arg(flag, value, req)
    if Argument.req == True:
        if Argument.flag not in sys.argv:
            print(errormsg)
            sys.exit("ERROR: required argument '{0}' is missing".format(Argument.flag))
    if Argument.value == True:
        try:
            test = sys.argv[sys.argv.index(Argument.flag)+1]
        except ValueError:
            if Argument.req == True:
                print(errormsg)
                sys.exit("ERROR: required argument '{0}' is missing".format(Argument.flag))
            elif Argument.req == False:
                return False
        except IndexError:
                print(errormsg)
                sys.exit("ERROR: argument '{0}' requires a value".format(Argument.flag))
        else:
            if Argument.value == True:
                Argument.value = sys.argv[sys.argv.index(Argument.flag)+1]
        
    if Argument.value == False:
        if Argument.flag in sys.argv:
            Argument.value = True
        else:
            Argument.value = False
    return Argument.value
#-----------------------------------------------------------------------------------#


def azimuthalAverage(image, center=None):
    #print('calculating rotationally avgeraged 1D PS')
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

    r = np.hypot(x - center[(0)], y - center[(1)])

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[(ind)]
    i_sorted = image.flat[(ind)]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[(0)]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin
    
    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[(rind[1:])] - csim[(rind[:-1])]

    radial_prof = tbin / nr

    return radial_prof

def get_PS(image_in):
    #print('calculating FFT and PS')
    mrcim = mrc_image(image_in)
    image = read(mrcim)
    #image = plt.imread(image_in)           # for reading from tiffs instead of mrcs

    # fourier transform of the image.

    F1 = fftpack.fft2(image)
     
    # make pretty - put low res in center
    F2 = fftpack.fftshift( F1 )
     
    # Calculate 2D power spectrum
    psd2D = np.abs( F2 )**2
     
    ## Calculate the azimuthally averaged 1D power spectrum
    psd1D = azimuthalAverage(psd2D)

    return (psd1D,psd2D,image)


def calc_fit(RAPS,PS2d,pxsize):
    samples = []
    for i in range(50,66):
        samples.append(int(round(calc_dist(PS2d,pxsize,i*0.1),0)))
    for i in range(28,42):
        samples.append(int(round(calc_dist(PS2d,pxsize,i*0.1),0)))
    
    values  = []
    for i in samples:
        values.append(np.log10(RAPS[(i)]))
    sams = np.asarray(samples)
    vals = np.asarray(values)
    
    fit = np.polyfit(sams, vals, 2)
    return(fit)

def calc_dist(ps2d,pxsize,resolution):
    """ given a resolution calculate its frequency
    on the RAPS
    """
    nyquist = 2*pxsize
    distance = ((np.shape(ps2d)[(0)]*.5)*nyquist)/resolution
    return(distance)

def make_plots(raps,ps2d,fact1,fact2,yint,pxsize,imagedata):
    fr1 = plt.subplot2grid((1,3), (0, 0),rowspan=2)

    fr1.axes.get_yaxis().set_visible(False)
    plt.plot(np.log10(raps))
    # calculate range for fit
    xmin = int(round(calc_dist(ps2d,pxsize,8.0),0))
    xmax = int(round(calc_dist(ps2d,pxsize,2.8),0))
    x = np.arange(xmin,xmax)
    eq = '{0}*x**2+{1}*x+{2}'.format(fact1,fact2,yint)
    #print('fit = {}'.format(eq))
    y = eval(eq)
    plt.plot(x,y,color='red',alpha=0.8)
    scatterpoints = [calc_dist(ps2d,pxsize,4.6),calc_dist(ps2d,pxsize,4.9)]
    vals = [(fact1*(scatterpoints[(0)]**2))+(fact2*scatterpoints[(0)])+yint,fact1*scatterpoints[(1)]**2+fact2*scatterpoints[(1)]+yint]
    plt.axvline(x=calc_dist(ps2d,pxsize,4.7),c='green', alpha=0.3)
    plt.axvline(x=calc_dist(ps2d,pxsize,4.9),c='green', alpha=0.3)
    
    fr1 = plt.subplot2grid((2,3), (0, 1),rowspan=2)
    fr1.axes.get_yaxis().set_visible(False)
    plt.xlim([calc_dist(ps2d,pxsize,6.0),calc_dist(ps2d,pxsize,4.0)])
    plt.plot(np.log10(raps))
    # calculate range for fit
    xmin = int(round(calc_dist(ps2d,pxsize,8.0),0))
    xmax = int(round(calc_dist(ps2d,pxsize,2.8),0))
    x = np.arange(xmin,xmax)
    eq = '{0}*x**2+{1}*x+{2}'.format(fact1,fact2,yint)
    #print('fit = {}'.format(eq))
    y = eval(eq)
    plt.plot(x,y,color='red',alpha=0.4)
    plt.axvline(x=calc_dist(ps2d,pxsize,4.6),c='green', alpha=0.3)
    plt.axvline(x=calc_dist(ps2d,pxsize,4.9),c='green', alpha=0.3)
    
    fr2 = plt.subplot2grid((2,3), (1, 2))

    fr2.axes.get_xaxis().set_visible(False)
    fr2.axes.get_yaxis().set_visible(False)
    plt.xlim(0,np.shape(ps2d)[1])
    plt.ylim(0,np.shape(ps2d)[0])
    plt.imshow(np.log10(ps2d),cmap='Greys_r')
    plt.plot((0.5*np.shape(ps2d)[0],0.5*np.shape(ps2d)[0]),(0.5*np.shape(ps2d)[0]+scatterpoints[0],0.5*np.shape(ps2d)[0]+scatterpoints[1]),color='orange')
    fr3 = plt.subplot2grid((2,3), (0, 2))
    fr3.axes.get_xaxis().set_visible(False)
    fr3.axes.get_yaxis().set_visible(False)
    plt.imshow(imagedata,cmap='Greys_r')
    
    return(plt)

def compare_vals(PS2d,raps,pxsize,fit):
    checkvals = []
    differences = []
    for i in range(470,490):
        checkvals.append(int(round(calc_dist(PS2d,pxsize,i*0.01),0)))
    for i in checkvals:
        calc = (fit[(0)]*i**2)+(fit[(1)]*i)+fit[(2)]
        actual = np.log10(raps[(i)])
        differences.append(actual - calc)
    return(np.mean(differences))

def do_it(image,pixelsize):
    imout = image.split('/')[-1]
    rotavg_PS,PS,imgdata = get_PS(image)
    fit = calc_fit(rotavg_PS,PS,pixelsize)
    plot = make_plots(rotavg_PS,PS,fit[(0)],fit[(1)],fit[(2)],pixelsize,imgdata)
    mean_dif = compare_vals(PS,rotavg_PS,pixelsize,fit)
    # calculate the expected values for 4.9 - 4.6 A
    plot.suptitle('score = {0}'.format(mean_dif))
    if dontwriteimgs == False:
        plot.savefig('fibril_quality_data/{0}_FQ.png'.format(imout.split('.')[0]))
    plot.close()
    return(mean_dif) 

### DO IT
#print("""
#-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=
#Fibril Quality Assessment v{0}
#-=-=-=--=-=-=-=-=-=-=-=-=-=-=-=""".format(vers))
## setuop


stack = make_arg('--stack',False,False)
imagesearch = make_arg('--i',True,True)
pixelsize = float(make_arg('--apix',True,True))
dontwriteimgs = make_arg('--noimages',False,False)

if dontwriteimgs == False:
    if os.path.isdir('fibril_quality_data') == False:
        os.system('mkdir fibril_quality_data')

if pixelsize > 2.0:
    sys.exit('ERROR: Pixel size must be less than 2.0 A/px to resolve cross-beta layerlines.\nERROR: Specified A/px == {0} '.format(pixelsize))
if imagesearch.split('.')[-1] in ['mrc','mrcs']:
    images = [imagesearch]
else:
    images = glob.glob('{0}*.mrc'.format(imagesearch))
imgs_dic = {}
means = {}
n=0
for i in images:
    means[i]=do_it(i,pixelsize)
    n+=1
    if n%10 == 0:
        sys.stdout.write('.')
        sys.stdout.flush()
# write output
output = open('FQ_find.log','a')
for i in means:
    output.write('{0}\t{1}\n'.format(i,means[i]))
output.close()

