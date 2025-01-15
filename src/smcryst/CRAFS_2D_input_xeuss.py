# -*- coding: utf-8 -*-
"""
Absorption correction not implemented correctlyT
"""

import pyFAI
from pyFAI.azimuthalIntegrator import AzimuthalIntegrator
import fabio
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as ticker
import argparse as ap



parser = ap.ArgumentParser(description=
                           "Perform data reduction on data from Xeuss 2.0C")

parser.add_argument("sample", help="sample file (.edf)")

parser.add_argument("mask", help="mask file (.edf)")
#parser.add_argument("poni", help="PONI configuration file (.poni)")
parser.add_argument("--empty", help="empty file (.edf)")
parser.add_argument("--background", help="background file (.edf)")
parser.add_argument("--opencl", help="Use OpenCL (requires pyopencl)",action="store_true")
parser.add_argument("--bkg_factor", help="background correction factor",type=float,default=float(1))
parser.add_argument("--eta", help="number of data points in eta",type=int,default=360)
parser.add_argument("--theta", help="number of data points in 2theta",type=int,default=2500)
parser.add_argument("--noabs", help="No absorption correction",action="store_true")
parser.add_argument("--autobkg", help="Automatic background scaling",action="store_true")
parser.add_argument("--do1D", help="output 1D",action="store_true")
parser.add_argument("--scale", help="Scale factor (used to help with fitting 1D)",type=float,default=float(1))
args = parser.parse_args()


file = fabio.open(args.sample)
mask = fabio.open(args.mask)

args.noabs = True

if args.noabs == False:
    def abscor (tth,intensity,trans):
        fabs = (np.log(trans)*((1/np.cos(np.radians(tth)))-1))/(np.exp(np.log(trans)*((1/np.cos(np.radians(tth)))-1))-1)
        cor_int = intensity/fabs
        return cor_int
else:
    def abscor (tth,intensity,trans):
        return intensity

if args.opencl==True:
    intmeth='opencl'
else:
    intmeth='csr'



geo = {               
      "detector": pyFAI.detectors.Detector(float(file.header['PSize_2']),float(file.header['PSize_1'])),
      "dist": float(file.header['SampleDistance']),
      "poni1": float(file.header['PSize_2'])*float(file.header['Center_2']),
      "poni2": float(file.header['PSize_1'])*float(file.header['Center_1']),
      "rot1": 0,
      "rot2": 0,
      "rot3": 0,
      "wavelength": float(file.header['WaveLength'])
}

poni = AzimuthalIntegrator(**geo)

transmission_sample = float(file.header['Transmission'])
print("Sample transmission:"+str(round(transmission_sample,4)))

#norm_sample = float(file.header['Intensity1'])
norm_sample = 1

data_2D = file.data/transmission_sample

if args.background != None:
    bkg = fabio.open(args.background)
    if np.allclose(float(bkg.header['Transmission']),1,atol=5e-3) == True:
        transmission_bkg = 1
    else:
        transmission_bkg = float(bkg.header['Transmission'])

    print("Background transmission:"+str(round(transmission_bkg,4)))

    #norm_bkg = float(bkg.header['Intensity1'])
    norm_bkg = 1
    
    bkg_2D = bkg.data/transmission_bkg
    
else:
    bkg_2D = 0
if args.empty != None:
    empty = fabio.open(args.empty)
    if np.allclose(float(empty.header['Transmission']),1,atol=5e-3) == True:
        transmission_empty = 1
    else:
        transmission_empty = float(empty.header['Transmission'])
    norm_empty = float(empty.header['Intensity1'])
    norm_empty = 1
    empty_2D = empty.data/transmission_empty
    if args.background != None:
        data_2I = (data_2D-empty_2D)-(bkg_2D-empty_2D)
    else:
        data_2I = data_2D-empty_2D
else:
    if args.background == None:
        data_2I = data_2D
    else:
        try:
            data_2I = data_2D - bkg_2D
        except ValueError:
            print("Incorrect background, ignoring")
            data_2I = data_2D

if args.do1D == True:
    data_1D = poni.integrate1d(file.data,
                               1000, correctSolidAngle=True,
                                      method=intmeth, radial_range=(None), azimuth_range=(None),
                                      unit="2th_deg", mask=mask.data, normalization_factor=norm_sample,
                                      metadata=None,error_model="poisson")
    q2,I2,sig2 = data_1D
    eta2 = 0
else:
    data_2I2 = poni.integrate2d(data_2I,args.theta,args.eta,mask=mask.data,unit="2th_deg")
    
    I2,q2,eta2 = data_2I2

if np.min(I2) <= 0:
    I2 -= np.min(I2)-1e-9

"""Polarization correction according to doi:10.1107/S0021889813014805"""
polfact = 0.5*(1+((np.cos(np.radians(q2)))**2))

I22 = np.multiply(I2,polfact)
I22 *= args.scale

output_head = np.insert(q2,0,float(file.header['WaveLength']))
if args.do1D == True:
    output_matrix = np.insert(I22,0,eta2)
    output_matrix = np.vstack((output_head,output_matrix))
else:
    output_matrix = np.insert(I22,0,eta2,axis=1)
    output_matrix = np.insert(output_matrix,0,output_head,axis=0)

np.savetxt(args.sample[:-3]+"csv", output_matrix,delimiter=',')