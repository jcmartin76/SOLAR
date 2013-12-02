# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

# ..ooOO0OOoo....ooOO0OOoo....ooOO0OOoo....ooOO0OOoo....ooOO0OOoo..
# Program wich calculates Energy aociated to the WLE in a solar flare
# December 2013
# JCMO & Milo-BC
# email: jcbuitragoc@unal.edu.co

# Import all pakages needed to run the code
import sunpy
import pyfits
import numpy
import matplotlib.pyplot as plt
from scipy import *
from pylab import *

# <codecell>

# Reading files
# We need to put the next three files into the 'input_files' folder.

# Defining names to the files
filename_data = 'input_files/CLEAN_INTEN_GONG_X_T.fits' # cube of intensigrams
filename_diff = 'input_files/CLEAN_INTEN_DIFF.fits' # cube of differences
filename_mask = 'input_files/mask_qs.fits' # mask to be applied

# Reading the files
hdu_gong = pyfits.open(filename_data)
hdu_diff = pyfits.open(filename_diff)
hdu_mask = pyfits.open(filename_mask)
 
# Each name is a numerical object, we devide data and header as:
datag = hdu_gong[0].data
hdrg  = hdu_gong[0].header

datad = hdu_diff[0].data
hdrd  = hdu_diff[0].header

datam = hdu_mask[0].data

# <codecell>


#LASP Irradiance of the day
lasp_irrad = 1357.50

fmask = 27

a = datag[fmask,:,:]

#Mask definition using the 70% contour
datar = numpy.ma.masked_where(datad[fmask,:,:] < 0.7*numpy.ma.maximum(datad[fmask,:,:]), a)


#Show masked image (opc.)

plt.imshow(datar, cmap=plt.cm.gray)
plt.show()


import os
import re 
dirList=os.listdir('RAW_DATA/')

total_Int_Sun = []
tt = []
dt = []

for fname in dirList:
    f = os.popen("sig_stats RAW_DATA/%s" %(fname),'r').readlines()#  > tmp_stats")

af = f[2].split()
mean_intensity = float(af[2])
af = f[3].split()
ntotal_pixels  =  float(af[0])

#Reading fits files to generate the temporal array
hdu_tmp = pyfits.open("RAW_DATA/%s" %(fname))
hdat  = hdu_tmp[0].data
hdat = rot90(hdat)
htmp  = hdu_tmp[0].header

qs_norm = numpy.mean(numpy.ma.masked_where(datam == 1.0, hdat))


str = "T"
seq = (htmp['DATE-OBS'], htmp['TIME-OBS']) # This is sequence of strings.
tt.append(str.join(seq))
    
str = " "
seq = (htmp['DATE-OBS'], htmp['TIME-OBS']) # This is sequence of strings.
dt.append(str.join(seq))

 
#Computation and storage of the total solar intensity
total_Int_Sun.append(mean_intensity*ntotal_pixels/qs_norm)


#------------------------------
#  Gong pixel in arcsec = 2.4
#  Gong pixel in Mm = 1.74
#------------------------------

#Image statitics
mask_stats = float(numpy.ma.count(datar))

print "Mask Image Statistics"  
print "Total of Pixel Values = %f" % (mask_stats)   
print "Total Area in px^2    = %f" % (numpy.sqrt(mask_stats))
print "Total Area in arcsec^2= %f" % (numpy.sqrt(mask_stats)*2.4)
print "Total Area in cm^2    = %g" % (numpy.sqrt(mask_stats)*1.74E16)
print " "
print "Total Solar Intensity = %g" % (numpy.mean(total_Int_Sun)/ntotal_pixels) 

pos_px_size = 1.0# (700.0E6 * 0.002)**2

#Computation of the irradiance as the Sum of all values inside the area.
irrad = []
ttp = []

# <codecell>

for i in range (0,len(datag[:,0,0])): 
    image = datag[i,:,:]
    I_wl = 2.0*pos_px_size*(numpy.sum(numpy.ma.masked_where(datad[fmask,:,:] < 0.7*numpy.ma.maximum(datad[fmask,:,:]), image))) #/numpy.mean(total_Int_Sun)
    irrad.append(I_wl*lasp_irrad/total_Int_Sun) 
    ttp.append(i)

yerr = 1.e-5   

# define our fitting function
from scipy.optimize import curve_fit

tarray = array(ttp)
iarray = array(irrad)

def fit_func(x, a0, a1, a2, a3, a4, a5, a6, a7):
    a1=fmask+1
    z = (x - a1) / a2
    stp = a5*(1+(2./numpy.pi)*numpy.arctan(a6*(x-a1)))
    y = a0 * numpy.exp(-z**2 / a2) + a3 + a4 * x  + stp + a7 * x**2
    return y


#Computation of the Luminosity
import ephem

u = ephem.Sun(dt[0])
eph  = u.earth_distance
L_WL = numpy.pi*eph*(149598000000.**2)*iarray #*lasp_irrad    #;*cos(!dtor*20.)

#Estimation of the error based on the preflare variation of the luminosity 
pf_max =  L_WL[0:20].max()
pf_min =  L_WL[0:20].min()

#Computation of the energy in ergs

bw_e = 0.0

# <codecell>

#Es necesario mirar la difinición de parameter y wfit

# <codecell>


# <codecell>


