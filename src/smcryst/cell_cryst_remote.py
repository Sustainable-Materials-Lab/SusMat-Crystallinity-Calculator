#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
Created on Sat Feb 28 2015 at 14:30:22

Python script for the generation of XRD graphs from TOPAS output and
generation of crystallinity data
"""
import os
import os.path
import copy
import argparse as ap
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import scipy.integrate as intg
from irods.session import iRODSSession
from irods.path import iRODSPath

env_file = os.getenv('IRODS_ENVIRONMENT_FILE', os.path.expanduser(
    '~/.irods/irods_environment.json'))

mpl.rcParams['font.family'] = 'sans-serif'

def cli():
    """Calculate the crystallinity of the cellulose sample and plot the data"""
    parser = ap.ArgumentParser(description='Plot XRD from topas academic')
    parser.add_argument("input", help="Path to iRODS collection containing TOPAS txt files")
    parser.add_argument(
        "--mono", help="Monochromated copper Kalpha used", action="store_true")
    parser.add_argument(
        "--png", help="png output (do not combine with svg)", action="store_true")
    parser.add_argument(
        "--svg", help="svg output (do not combine with png)", action="store_true")
    parser.add_argument("--onecol", help="single column output",
                        action="store_true")
    parser.add_argument(
        "--clip", help="Clip 2theta range to xmin-xmax", action="store_true")
    parser.add_argument(
        "--celtype", help="type of cellulose(use latex math notation)", default=r"Cellulose I$_\beta$")
    parser.add_argument(
        "--cel2type", help="type of cellulose 2 (use latex math notation)", default=r"Cellulose II")
    parser.add_argument(
        "--ymax", help="Y max offset from data (percent)", type=float, default=1)
    parser.add_argument(
        "--xmax", help="Max 2theta used for clipping", type=float, default=55)
    parser.add_argument(
        "--xmin", help="Min 2theta used for clipping", type=float, default=5)
    parser.add_argument(
        "--legloc", help="legend location eg. upper left", type=str, default="upper left")
    parser.add_argument(
        "--exposure", help="convert CPS input to counts (exp in seconds)", type=float, default=1)
    parser.add_argument(
        "--cel2", help="cellulose II phase present", action="store_true")
    parser.add_argument(
        "--linsub", help="Subtract linear portion of amorphous scattering", action="store_true")
    parser.add_argument(
        "--xye", help="input file had errors included", action="store_true")
    parser.add_argument(
        "--keepbkg", help="Do not subtract the background", action="store_true")
    parser.add_argument(
        "--peaks", help="amorphous phase consists of xo_Is peaks", action="store_true")
    parser.add_argument("--raw", help="Also plot raw data", action="store_true")
    parser.add_argument(
        "--bckamorph", help="Background models amorphous content", action="store_true")
    parser.add_argument("--theta", help="plot in 2theta", action="store_true")
    parser.add_argument(
        "--absolute", help="use absolute intensity", action="store_true")
    parser.add_argument(
        "--rawcryst", help="use raw data for crystallinity (rather than profile), use total cellulose for cryst. phase", action="store_true")
    parser.add_argument(
        "--FOLorentz", help="No Lorentz correction", action="store_true")

    args = parser.parse_args()  # map arguments to args

    if args.mono is True:
        lambda_1 = 1.54056
    else:
        lambda_1 = 1.54189


    width = 17.78

    if args.onecol is True:
        width = width/2
    winch = width/2.54  # pyplot is american so uses inches

    hinch = winch/1.75  # ISO aspect ratio

    if width <= 9:  # Change font sizes based on image size
        labsize = 6
        legsize = 4
    elif width > 9 and width <= 15:
        labsize = 8
        legsize = 6
    else:
        labsize = 10
        legsize = 6

    # Ask for filetype
    if args.png is True:
        filetp = 'png'
    elif args.svg is True:
        filetp = 'svg'
    else:
        filetp = 'pdf'

    dpipng = 600

    with iRODSSession(irods_env_file=env_file) as session:
        coll = session.collections.get(iRODSPath(args.input))
        for obj in coll.data_objects:
            if obj.name.endswith(".txt"):
                with obj.open("r") as infile:
                    rawtxt = infile.read().decode('utf-8').splitlines()
                    data = np.genfromtxt(rawtxt, delimiter=',', skip_header=2,
                                        invalid_raise=False, unpack=False)

                    #Detect the presence of extra phases

                    if "Sqrt(y)" in rawtxt[0]:
                        args.sqrt = True
                        data2 = np.square(data)
                        data2 = np.delete(data2, 0, 1)
                        data = np.insert(data2, 0, data[:, 0], 1)


                    if "Alpha i-PP" in rawtxt[1]:
                        args.iPP = True
                    else:
                        args.iPP = False

                    if "Gamma i-PP" in rawtxt[1]:
                        args.giPP = True
                    else:
                        args.giPP = False

                    if "PCL" in rawtxt[1]:
                        args.PCL = True
                    else:
                        args.PCL = False

                    if "Cellulose II" in rawtxt[1]:
                        args.cel2 = True
                    elif args.cel2 is True:
                        args.cel2 = True
                    else:
                        args.cel2 = False

                    if "Amorphous,Background" in rawtxt[1]:
                        args.cel2amorph = True
                    else:
                        args.cel2amorph = False
                    if "Bkg" in rawtxt[1]:
                        args.background = True
                    else:
                        args.background = False
                    if "Background,Amorphous" in rawtxt[1]:
                        args.expback = True
                    else:
                        args.expback = False

                    data = data.T

                    ttheta = []

                    #Import Raw Data
                    if args.peaks is True:
                        if args.xye is True:
                            if args.cel2 is True:
                                ttheta, raw, _, prf, diff, cel1, cel2, = data

                                amorph = prf - cel1 - cel2
                            elif args.iPP is True:
                                if args.giPP is True:
                                    ttheta, raw, _, prf, diff, cel1, iPP, gipp = data

                                    amorph = prf - cel1 - gipp - iPP
                                else:
                                    ttheta, raw, _, prf, diff, cel1, iPP = data

                                    amorph = prf - cel1 - iPP
                            elif args.PCL is True:
                                ttheta, raw, _, prf, diff, cel1, PCL = data

                                amorph = prf - cel1 - PCL
                            else:
                                ttheta, raw, _, prf, diff, cel1 = data

                                amorph = prf - cel1
                        else:
                            if args.cel2 is True:
                                ttheta, raw, prf, diff, cel1, cel2 = data

                                amorph = prf - cel1 - cel2
                            elif args.iPP is True:
                                if args.giPP is True:
                                    ttheta, raw, prf, diff, cel1, iPP, gipp = data

                                    amorph = prf - cel1 - gipp - iPP
                                else:
                                    ttheta, raw, prf, diff, cel1, iPP = data

                                    amorph = prf - cel1 - iPP
                            elif args.PCL is True:
                                ttheta, raw, prf, diff, cel1, PCL = data

                                amorph = prf - cel1 - PCL
                            else:
                                ttheta, raw, prf, diff, cel1 = data

                                amorph = prf - cel1
                    elif args.background is True:
                        if args.bckamorph is True:
                            if args.xye is True:
                                if args.cel2 is True:
                                    ttheta, raw, _, bck, prf, diff, cel1, cel2 = data

                                    cel2 -= bck
                                elif args.iPP is True:
                                    if args.giPP is True:
                                        ttheta, raw, _, bck, prf, diff, cel1, iPP, gipp = data

                                        iPP -= bck
                                        gipp -= bck
                                    else:
                                        ttheta, raw, _, bck, prf, diff, cel1, iPP = data

                                        iPP -= bck
                                elif args.PCL is True:
                                    ttheta, raw, _, bck, prf, diff, cel1, PCL = data

                                    PCL -= bck
                                else:
                                    ttheta, raw, _, bck, prf, diff, cel1 = data

                            else:
                                if args.cel2 is True:
                                    ttheta, raw, bck, prf, diff, cel1, cel2 = data

                                    cel2 -= bck
                                elif args.iPP is True:
                                    if args.giPP is True:
                                        ttheta, raw, bck, prf, diff, cel1, iPP, gipp = data

                                        iPP -= bck
                                        gipp -= bck
                                    else:
                                        ttheta, raw, bck, prf, diff, cel1, iPP = data

                                        iPP -= bck
                                elif args.PCL is True:
                                    ttheta, raw, bck, prf, diff, cel1, PCL = data

                                    PCL -= bck
                                else:
                                    ttheta, raw, bck, prf, diff, cel1 = data

                            cel1 -= bck
                            amorph = copy.copy(bck)
                            bck *= 0
                        else:
                            if args.xye is True:
                                if args.cel2 is True:
                                    ttheta, raw, _, bck, prf, diff, cel1, cel2, amorph = data

                                    prf -= bck
                                    raw -= bck
                                    amorph -= bck
                                    cel1 -= bck
                                    cel2 -= bck
                                elif args.iPP is True:
                                    if args.giPP is True:
                                        ttheta, raw, _, bck, prf, diff, cel1, iPP, gipp, amorph = data

                                        prf -= bck
                                        amorph -= bck
                                        cel1 -= bck
                                        iPP -= bck
                                        gipp -= bck
                                        raw -= bck
                                    else:
                                        ttheta, raw, _, bck, prf, diff, cel1, iPP, amorph = data

                                        prf -= bck
                                        amorph -= bck
                                        cel1 -= bck
                                        iPP -= bck

                                        raw -= bck
                                elif args.PCL is True:
                                    ttheta, raw, _, bck, prf, diff, cel1, PCL, amorph = data

                                    prf -= bck
                                    amorph -= bck
                                    cel1 -= bck
                                    PCL -= bck

                                    raw -= bck
                                else:
                                    ttheta, raw, _, bck, prf, diff, cel1, amorph = data

                                    prf -= bck
                                    amorph -= bck
                                    cel1 -= bck
                                    raw -= bck
                            else:
                                if args.cel2 is True:
                                    ttheta, raw, bck, prf, diff, cel1, cel2, amorph = data

                                    prf -= bck
                                    amorph -= bck
                                    cel1 -= bck
                                    cel2 -= bck
                                    amorph -= bck
                                    raw -= bck
                                elif args.iPP is True:
                                    if args.giPP is True:
                                        ttheta, raw, bck, prf, diff, cel1, iPP, gipp, amorph = data

                                        prf -= bck
                                        amorph -= bck
                                        cel1 -= bck
                                        iPP -= bck
                                        gipp -= bck
                                        raw -= bck
                                    else:
                                        ttheta, raw, bck, prf, diff, cel1, iPP, amorph = data

                                        prf -= bck
                                        amorph -= bck
                                        cel1 -= bck
                                        iPP -= bck
                                        amorph -= bck
                                        raw -= bck
                                elif args.PCL is True:
                                    ttheta, raw, bck, prf, diff, cel1, PCL, amorph = data

                                    prf -= bck
                                    amorph -= bck
                                    cel1 -= bck
                                    PCL -= bck
                                    amorph -= bck
                                    raw -= bck
                                else:
                                    ttheta, raw, bck, prf, diff, cel1, amorph = data

                                    prf -= bck
                                    amorph -= bck
                                    cel1 -= bck
                                    amorph -= bck
                                    raw -= bck
                    elif args.expback is True:
                        if args.xye is True:
                            if args.cel2 is True:
                                ttheta, raw, _, prf, diff, cel1, cel2, bck, amorph = data

                                prf -= bck
                                raw -= bck
                                amorph -= bck
                                cel1 -= bck
                                cel2 -= bck
                            elif args.iPP is True:
                                if args.giPP is True:
                                    ttheta, raw, _, prf, diff, cel1, iPP, gipp, bck, amorph = data

                                    prf -= bck
                                    amorph -= bck
                                    cel1 -= bck
                                    iPP -= bck
                                    gipp -= bck
                                    raw -= bck
                                else:
                                    ttheta, raw, _, prf, diff, cel1, iPP, bck, amorph = data

                                    prf -= bck
                                    amorph -= bck
                                    cel1 -= bck
                                    iPP -= bck
                                    raw -= bck
                            elif args.PCL is True:
                                ttheta, raw, _, prf, diff, cel1, PCL, bck, amorph = data

                                prf -= bck
                                amorph -= bck
                                cel1 -= bck
                                PCL -= bck
                                raw -= bck
                            else:
                                ttheta, raw, _, prf, diff, cel1, bck, amorph = data

                                prf -= bck
                                amorph -= bck
                                cel1 -= bck
                                raw -= bck
                        else:
                            if args.cel2 is True:
                                ttheta, raw, prf, diff, cel1, cel2, bck, amorph = data

                                prf -= bck
                                amorph -= bck
                                cel1 -= bck
                                cel2 -= bck
                                amorph -= bck
                                raw -= bck
                            elif args.iPP is True:
                                if args.giPP is True:
                                    ttheta, raw, prf, diff, cel1, iPP, gipp, bck, amorph = data

                                    prf -= bck
                                    amorph -= bck
                                    cel1 -= bck
                                    iPP -= bck
                                    gipp -= bck
                                    raw -= bck
                                else:
                                    ttheta, raw, prf, diff, cel1, iPP, bck, amorph = data

                                    prf -= bck
                                    amorph -= bck
                                    cel1 -= bck
                                    iPP -= bck
                                    amorph -= bck
                                    raw -= bck
                            elif args.PCL is True:
                                ttheta, raw, prf, diff, cel1, PCL, bck, amorph = data

                                prf -= bck
                                amorph -= bck
                                cel1 -= bck
                                PCL -= bck
                                amorph -= bck
                                raw -= bck
                            else:
                                ttheta, raw, prf, diff, cel1, bck, amorph = data

                                prf -= bck
                                amorph -= bck
                                cel1 -= bck
                                amorph -= bck
                                raw -= bck
                    else:
                        if args.xye is True:
                            if args.cel2 is True:
                                if args.cel2amorph is True:
                                    ttheta, raw, _, prf, diff, cel1, cel2, amorph = data

                                else:
                                    ttheta, raw, _, prf, diff, cel1, cel2, amorph = data

                            elif args.iPP is True:
                                if args.giPP is True:
                                    if args.cel2amorph is True:
                                        ttheta, raw, _, prf, diff, cel1, iPP, amorph, gipp = data

                                    else:
                                        ttheta, raw, _, prf, diff, cel1, iPP, amorph, gipp = data

                                else:
                                    if args.cel2amorph is True:
                                        ttheta, raw, _, prf, diff, cel1, iPP, amorph = data

                                    else:
                                        ttheta, raw, _, prf, diff, cel1, iPP, amorph = data

                            elif args.PCL is True:
                                if args.cel2amorph is True:
                                    ttheta, raw, _, prf, diff, cel1, PCL, amorph = data

                                else:
                                    ttheta, raw, _, prf, diff, cel1, PCL, amorph = data

                            else:
                                if args.cel2amorph is True:
                                    ttheta, raw, _, prf, diff, cel1, amorph = data
                                    # Raw data
                                else:
                                    ttheta, raw, _, prf, diff, cel1, amorph = data

                        else:
                            if args.cel2 is True:
                                if args.cel2amorph is True:
                                    ttheta, raw, prf, diff, cel1, cel2, amorph = data

                                else:
                                    ttheta, raw, prf, diff, cel1, cel2, amorph = data

                            elif args.iPP is True:
                                if args.giPP is True:
                                    if args.cel2amorph is True:
                                        ttheta, raw, prf, diff, cel1, iPP, amorph, gipp = data

                                    else:
                                        ttheta, raw, prf, diff, cel1, iPP, amorph, gipp = data

                                else:
                                    if args.cel2amorph is True:
                                        ttheta, raw, prf, diff, cel1, iPP, amorph = data

                                    else:
                                        ttheta, raw, prf, diff, cel1, iPP, amorph = data

                            elif args.PCL is True:
                                if args.cel2amorph is True:
                                    ttheta, raw, prf, diff, cel1, PCL, amorph = data

                                else:
                                    ttheta, raw, prf, diff, cel1, PCL, amorph = data

                            else:
                                if args.cel2amorph is True:
                                    ttheta, raw, prf, diff, cel1, amorph = data
                                    # Raw data
                                else:
                                    ttheta, raw, prf, diff, cel1, amorph = data

                    raw *= args.exposure
                    prf *= args.exposure
                    diff *= args.exposure
                    cel1 *= args.exposure
                    amorph *= args.exposure
                    if args.cel2 is True:
                        cel2 *= args.exposure
                    if args.iPP is True:
                        iPP *= args.exposure
                    if args.giPP is True:
                        gipp *= args.exposure
                    if args.PCL is True:
                        PCL *= args.exposure

                    if args.clip is True:
                        startval = np.where(ttheta > args.xmin)[0][0]
                        endval = np.where(ttheta > args.xmax)[0][0]
                        ttheta = ttheta[startval:endval]
                        raw = raw[startval:endval]
                        prf = prf[startval:endval]
                        diff = diff[startval:endval]
                        bck = bck[startval:endval]
                        cel1 = cel1[startval:endval]
                        amorph = amorph[startval:endval]
                        if args.cel2 is True:
                            cel2 = cel2[startval:endval]
                        if args.iPP is True:
                            iPP = iPP[startval:endval]
                        if args.giPP is True:
                            gipp = gipp[startval:endval]
                        if args.PCL is True:
                            PCL = PCL[startval:endval]

                    if args.linsub is True:  # Remove the linear portion of the "amorphous" diffraction and add to the background
                        minamorph = min(amorph)
                        amorph -= minamorph
                        mincel1 = min(cel1)
                        cel1 -= mincel1
                        if args.cel2 is True:
                            mincel2 = min(cel2)
                            cel2 -= mincel2
                            prf -= mincel2
                            raw -= mincel2
                        if args.iPP is True:
                            miniPP = min(iPP)
                            iPP -= miniPP
                            prf -= miniPP
                            raw -= miniPP
                        if args.giPP is True:
                            mingipp = min(gipp)
                            gipp -= mingipp
                            prf -= mingipp
                            raw -= mingipp
                        if args.PCL is True:
                            minpcl = min(PCL)
                            PCL -= minpcl
                            prf -= minpcl
                            raw -= minpcl
                        prf = prf - minamorph - mincel1
                        raw = raw - minamorph - mincel1
                    # remove background from cel1 and amorph
                    # For plotting add background first

                    if args.keepbkg is True:
                        raw += bck
                        cel1 += bck
                        prf += bck
                        amorph += bck
                        if args.cel2 is True:
                            cel2 += bck
                        if args.iPP is True:
                            iPP += bck
                        if args.giPP is True:
                            gipp += bck
                        if args.PCL is True:
                            PCL += bck

                    tot = prf

                    if args.rawcryst is True:
                        tot = raw
                        cryst = prf - amorph

                    # Copper Kalpha1 1.54056 Angstroms Kalpha=1.5418 see input black above
                    lambda_2 = lambda_1/10  # in nm
                    # s = ((2*np.sin((ttheta/2)*(np.pi/180)))/lambda_2) #numpy uses radians
                    s = ((4*np.pi*np.sin((ttheta/2)*(np.pi/180)))/lambda_2)  # is Q

                    if args.FOLorentz is True:
                        #iraw1 = intg.cumulative_trapezoid((raw), s, initial=0)  # Integral curve
                        #iraw2 = intg.trapezoid((raw), s)  # Integral value
                        icel12 = intg.trapezoid(cel1, s)
                        iamorph1 = intg.cumulative_trapezoid(
                            (amorph), s, initial=0)  # Integral curve
                        iamorph2 = intg.trapezoid((amorph), s)  # Integral value
                        #iprf1 = intg.cumulative_trapezoid((prf), s, initial=0)  # Integral curve
                        #iprf2 = intg.trapezoid((prf), s)  # Integral value
                        itot1 = intg.cumulative_trapezoid((tot), s, initial=0)  # Integral curve
                        itot2 = intg.trapezoid((tot), s)  # Integral value
                        if args.rawcryst is True:
                            icryst = intg.trapezoid((cryst), s)
                            icryst1 = intg.cumulative_trapezoid((cryst), s, initial=0)
                        else:
                            icryst = itot2-iamorph2
                            icryst1 = itot1 - iamorph1
                    else:
                        #iraw1 = intg.cumulative_trapezoid(
                        #    (raw*(s**2)), s, initial=0)  # Integral curve
                        #iraw2 = intg.trapezoid((raw*(s**2)), s)  # Integral value
                        icel12 = intg.trapezoid(cel1*(s**2), s)
                        iamorph1 = intg.cumulative_trapezoid(
                            (amorph*(s**2)), s, initial=0)  # Integral curve
                        iamorph2 = intg.trapezoid((amorph*(s**2)), s)  # Integral value
                        #iprf1 = intg.cumulative_trapezoid(
                        #   (prf*(s**2)), s, initial=0)  # Integral curve
                        #iprf2 = intg.trapezoid((prf*(s**2)), s)  # Integral value
                        itot1 = intg.cumulative_trapezoid(
                            (tot*(s**2)), s, initial=0)  # Integral curve
                        itot2 = intg.trapezoid((tot*(s**2)), s)  # Integral value
                        if args.rawcryst is True:
                            icryst = intg.trapezoid((cryst*(s**2)), s)
                            icryst1 = intg.cumulative_trapezoid((cryst*(s**2)), s, initial=0)
                        else:
                            icryst = intg.trapezoid(((tot-amorph)*(s**2)), s)
                            icryst1 = intg.cumulative_trapezoid(
                                ((tot-amorph)*(s**2)), s, initial=0)
                    chi_c = icryst/itot2  # Crystallinity index

                    fcel1 = icel12/icryst

                    if args.cel2 is True:
                        icel22 = intg.trapezoid(cel2*(s**2), s)
                        fcel2 = icel22/icryst
                    if args.iPP is True:
                        iiPP2 = intg.trapezoid(iPP*(s**2), s)
                        fiPP = iiPP2/icryst
                    if args.giPP is True:
                        igipp2 = intg.trapezoid(gipp*(s**2), s)
                        fgipp = igipp2/icryst
                    if args.PCL is True:
                        iPCL2 = intg.trapezoid(PCL*(s**2), s)
                        fPCL = iPCL2/icryst

                    # create figure and axes
                    fig, ax = plt.subplots(figsize=(winch, hinch), layout='constrained') #pylint: disable=unused-variable

                    # Cellulose 1
                    if args.theta is True:
                        ax1 = ax.scatter(ttheta, raw, s=0.5, color='#9D9D92')
                        ax11, = ax.plot(ttheta, (cel1), color='b', linewidth=.5)
                        if args.cel2 is True:
                            ax11a, = ax.plot(ttheta, cel2, color='g', linewidth=.5)
                        if args.iPP is True:
                            ax11a, = ax.plot(ttheta, iPP, color='g', linewidth=.5)
                        if args.giPP is True:
                            ax12a, = ax.plot(ttheta, gipp, color='orange', linewidth=.5)
                        if args.PCL is True:
                            ax11a, = ax.plot(ttheta, PCL, color='g', linewidth=.5)
                        # Amorph cellulose
                        ax12, = ax.plot(ttheta, (amorph), color='r', linewidth=.5)
                        # Background
                        # Profile
                        # comma required after ax2 so that legend works
                        ax2, = ax.plot(ttheta, prf)
                        plt.setp(ax2, linewidth=0.5, color='#000000')  # set colour black

                        # Plot on right axis
                        twinax = ax.twinx()
                        ax4, = twinax.plot(ttheta, itot1)
                        plt.setp(ax4, linewidth=0.5, linestyle='--', color='#000000')
                        ax5, = twinax.plot(ttheta, icryst1)
                        plt.setp(ax5, linewidth=0.5, linestyle='--', color='r')

                        # Dummy item for chi_c in legend
                        dummy = Rectangle((0, 0), 1, 1, fc="w", fill=False,
                                        edgecolor='none', linewidth=0)

                        # tick-marks
                        ax.tick_params(axis='both', which='major',
                                    labelsize=labsize)  # Axis number fontsize
                        ax.tick_params(axis='both', which='major',
                                    labelsize=labsize)  # Axis number fontsize
                        twinax.tick_params(axis='both', which='major', labelsize=labsize)
                        # Switch off top ticks, make bottom 'out'
                        ax.get_xaxis().set_tick_params(which='both', direction='out', top=False)
                        # Switch off right ticks on first set of axes, make bottom 'out'
                        ax.get_yaxis().set_tick_params(which='both', direction='out', right=False)
                        ax.ticklabel_format(style='sci', axis='y',
                                            scilimits=(0, 0), useMathText=True)
                        ax.yaxis.get_offset_text().set_size(labsize)
                        twinax.get_yaxis().set_tick_params(which='both', direction='out', right=True)
                        twinax.ticklabel_format(style='sci', axis='y',
                                                scilimits=(0, 0), useMathText=True)
                        twinax.yaxis.get_offset_text().set_size(labsize)
                        # Axes
                        plt.xlim(ttheta[0], ttheta[-1])  # X range
                        # y range based on max of data
                        ax.set_ylim(0, ((args.ymax/100)+1)*max(*prf, *raw))
                        twinax.set_ylim(0, (max(itot1)))
                        # Labels
                        ax.set_ylabel(r'$I$ /cm$^{-1}$', size=labsize)
                        # Can use LaTeX in labels
                        ax.set_xlabel(r'$2\theta$ /$^{\circ}$', size=labsize)
                        if args.FOLorentz is True:
                            twinax.set_ylabel(r'$\int I\left(q\right) dq$', size=labsize)
                        else:
                            twinax.set_ylabel(r'$\int I\left(q\right)q^2 dq$', size=labsize)
                    else:
                        ax1 = ax.scatter(s, raw, s=0.5, color='#9D9D92')
                        ax11, = ax.plot(s, (cel1), color='b', linewidth=.5)
                        if args.cel2 is True:
                            ax11a, = ax.plot(s, cel2, color='g', linewidth=.5)
                        if args.iPP is True:
                            ax11a, = ax.plot(s, iPP, color='g', linewidth=.5)
                        if args.giPP is True:
                            ax12a, = ax.plot(s, gipp, color='orange', linewidth=.5)
                        if args.PCL is True:
                            ax11a, = ax.plot(s, PCL, color='g', linewidth=.5)
                        # Amorph cellulose
                        ax12, = ax.plot(s, (amorph), color='r', linewidth=.5)
                        # Background
                        # Profile
                        ax2, = ax.plot(s, prf)  # comma required after ax2 so that legend works
                        plt.setp(ax2, linewidth=0.5, color='#000000')  # set colour black

                        # Plot on right axis
                        twinax = ax.twinx()
                        ax4, = twinax.plot(s, itot1)
                        plt.setp(ax4, linewidth=0.5, linestyle='--', color='#000000')
                        ax5, = twinax.plot(s, icryst1)
                        plt.setp(ax5, linewidth=0.5, linestyle='--', color='r')

                        # Dummy item for chi_c in legend
                        dummy = Rectangle((0, 0), 1, 1, fc="w", fill=False,
                                        edgecolor='none', linewidth=0)

                        # tick-marks
                        ax.tick_params(axis='both', which='major',
                                    labelsize=labsize)  # Axis number fontsize
                        ax.tick_params(axis='both', which='major',
                                    labelsize=labsize)  # Axis number fontsize
                        twinax.tick_params(axis='both', which='major', labelsize=labsize)
                        # Switch off top ticks, make bottom 'out'
                        ax.get_xaxis().set_tick_params(which='both', direction='out', top=False)
                        # Switch off right ticks on first set of axes, make bottom 'out'
                        ax.get_yaxis().set_tick_params(which='both', direction='out', right=False)
                        ax.ticklabel_format(style='sci', axis='y',
                                            scilimits=(0, 0), useMathText=True)
                        ax.yaxis.get_offset_text().set_size(labsize)
                        twinax.get_yaxis().set_tick_params(which='both', direction='out', right=True)
                        twinax.ticklabel_format(style='sci', axis='y',
                                                scilimits=(0, 0), useMathText=True)
                        twinax.yaxis.get_offset_text().set_size(labsize)
                        # Axes
                        plt.xlim(s[0], s[-1])  # X range
                        # y range based on max of data
                        ax.set_ylim(0, ((args.ymax/100)+1)*max(*prf, *raw))
                        twinax.set_ylim(0, (max(itot1)))
                        # Labels
                        if args.absolute is True:
                            ax.set_ylabel(r'$I$ /cm$^{-1}$', size=labsize)
                        else:
                            ax.set_ylabel(r'$I$ /Arb.', size=labsize)
                        # Can use LaTeX in labels
                        ax.set_xlabel(r'$q=4\pi\sin\theta/\lambda$ /nm$^{-1}$', size=labsize)
                        if args.FOLorentz is True:
                            twinax.set_ylabel(r'$\int I\left(q\right) dq$', size=labsize)
                        else:
                            twinax.set_ylabel(r'$\int I\left(q\right)q^2 dq$', size=labsize)


                    # Legend
                    #irawstr = str(round(float(iraw2), 2))  # integral to 2dp for label
                    iamorphstr = str(round(float(icryst), 2))
                    itotstr = str(round(float(itot2), 2))
                    chicstr = str(round(float(chi_c), 2))
                    cel1str = str(round(float(100*fcel1), 2))
                    if args.cel2 is True:
                        cel2str = str(round(float(100*fcel2), 2))
                    if args.iPP is True:
                        iPPstr = str(round(float(100*fiPP), 2))
                    if args.giPP is True:
                        gippstr = str(round(float(100*fgipp), 2))
                    if args.PCL is True:
                        PCLstr = str(round(float(100*fPCL), 2))

                    alpha = "\u03B1"
                    beta = "\u03B2"
                    gamma = "\u03B3"
                    chi = "\u03C7"


                    if args.cel2 is True:
                        # pylint: disable=used-before-assignment
                        plt.legend([ax1, ax11, ax11a, ax12, ax2, ax4, ax5, dummy],
                                ["Raw Data", f"{args.celtype}: {cel1str}%", f"{args.cel2type}: {cel2str}%", "Amorphous", "Fitted Profile", f"Profile Integral: {itotstr}",
                                    f"Crystalline Integral: {iamorphstr}", f"{chi}$_c$ = {chicstr}"],
                                loc=args.legloc, frameon=False, prop={'size': legsize})
                    if args.iPP is True:
                        if args.giPP is True:
                            # pylint: disable=used-before-assignment
                            plt.legend([ax1, ax11, ax11a, ax12a, ax12, ax2, ax4, ax5, dummy],
                                    ["Raw Data", f"Cellulose I{beta}: {cel1str}%", f"{alpha}-iPP: {iPPstr}%", f"{gamma}-iPP: {gippstr}%", "Amorphous", "Fitted Profile", f"Profile Integral: {itotstr}",
                                        f"Crystalline Integral: {iamorphstr}", f"{chi}$_c$ = {chicstr}"],
                                    loc=args.legloc, frameon=False, prop={'size': legsize})
                        else:
                            plt.legend([ax1, ax11, ax11a, ax12, ax2, ax4, ax5, dummy],
                                    ["Raw Data", f"Cellulose I{beta}: {cel1str}%", f"{alpha}-iPP: {iPPstr}%", "Amorphous", "Fitted Profile", f"Profile Integral: {itotstr}",
                                        f"Crystalline Integral: {iamorphstr}", f"{chi}$_c$ = {chicstr}"],
                                    loc=args.legloc, frameon=False, prop={'size': legsize})
                    if args.PCL is True:
                        plt.legend([ax1, ax11, ax11a, ax12, ax2, ax4, ax5, dummy],
                                ["Raw Data", f"Cellulose I{beta}: {cel1str}%", f"PCL: {PCLstr}%", "Amorphous", "Fitted Profile", f"Profile Integral: {itotstr}",
                                    f"Crystalline Integral: {iamorphstr}", f"{chi}$_c$ = {chicstr}"],
                                loc=args.legloc, frameon=False, prop={'size': legsize})
                    if (args.cel2 is False and args.iPP is False and args.PCL is False):
                        plt.legend([ax1, ax11, ax12, ax2, ax4, ax5, dummy],
                                ["Raw Data", args.celtype, "Amorphous", "Fitted Profile", f"Profile Integral: {itotstr}",
                                    f"Crystalline Integral: {iamorphstr}", f"{chi}$_c$ = {chicstr}"],
                                loc=args.legloc, frameon=False, prop={'size': legsize})


                    # Decide output
                    obj2 = session.data_objects.create(iRODSPath(args.input,obj.name[:-3]+filetp))
                    with obj2.open('r+') as outfile:
                        if filetp in ('png', 'PNG'):
                            plt.savefig(outfile, bbox_inches='tight', dpi=dpipng)
                        else:
                            plt.savefig(outfile, bbox_inches='tight')

                    if args.raw is True:
                        if args.theta is True:

                            # create figure and axes
                            fig, ax = plt.subplots(figsize=(hinch, winch), layout='constrained')

                            # Cellulose 1

                            ax1 = ax.plot(ttheta, raw, linewidth=0.5, color='k')
                            # Background

                            # tick-marks
                            ax.tick_params(axis='both', which='major',
                                        labelsize=labsize)  # Axis number fontsize
                            ax.tick_params(axis='both', which='major',
                                        labelsize=labsize)  # Axis number fontsize
                            # Switch off top ticks, make bottom 'out'
                            ax.get_xaxis().set_tick_params(which='both', direction='out', top=False)
                            # Switch off right ticks on first set of axes, make bottom 'out'
                            ax.get_yaxis().set_tick_params(which='both', direction='out', right=False)
                            ax.ticklabel_format(style='sci', axis='y',
                                                scilimits=(0, 0), useMathText=True)
                            ax.yaxis.get_offset_text().set_size(labsize)
                            # Axes
                            plt.xlim(ttheta[0], ttheta[-1])  # X range
                            # y range based on max of data
                            ax.set_ylim(0, ((args.ymax/100)+1)*max(*prf, *raw))
                            # Labels
                            ax.set_ylabel(r'$I$ /counts', size=labsize)
                            # Can use LaTeX in labels
                            ax.set_xlabel(r'$2\theta$ /$^{\circ}$', size=labsize)

                            # Legend

                            # A4 = 11.7,8.3 and single col = 3.25,1.5
                            plt.gcf().set_size_inches([winch, hinch])

                            # Decide output
                            obj2 = session.data_objects.create(iRODSPath(args.input,obj.name[:-4]+"_raw."+filetp))
                            with obj2.open('r+') as outfile:
                                if filetp in ('png', 'PNG'):
                                    plt.savefig(outfile, bbox_inches='tight', dpi=dpipng)
                                else:
                                    plt.savefig(outfile, bbox_inches='tight')
                        else:

                            # create figure and axes
                            fig, ax = plt.subplots(figsize=(hinch, winch), layout='constrained')
                            # Cellulose 1

                            ax1 = ax.plot(s, raw, linewidth=0.5, color='k')
                            # Background

                            # tick-marks
                            ax.tick_params(axis='both', which='major',
                                        labelsize=labsize)  # Axis number fontsize
                            ax.tick_params(axis='both', which='major',
                                        labelsize=labsize)  # Axis number fontsize
                            # Switch off top ticks, make bottom 'out'
                            ax.get_xaxis().set_tick_params(which='both', direction='out', top=False)
                            # Switch off right ticks on first set of axes, make bottom 'out'
                            ax.get_yaxis().set_tick_params(which='both', direction='out', right=False)
                            ax.ticklabel_format(style='sci', axis='y',
                                                scilimits=(0, 0), useMathText=True)
                            ax.yaxis.get_offset_text().set_size(labsize)
                            # Axes
                            plt.xlim(s[0], s[-1])  # X range
                            # y range based on max of data
                            ax.set_ylim(0, ((args.ymax/100)+1)*max(*prf, *raw))
                            # Labels
                            ax.set_ylabel(r'$I$ /counts', size=labsize)
                            # Can use LaTeX in labels
                            ax.set_xlabel(r'$q=4\pi\sin\theta/\lambda$ /nm$^{-1}$', size=labsize)

                            # Legend

                            # Decide output
                            obj2 = session.data_objects.create(iRODSPath(args.input,obj.name[:-4]+"_raw."+filetp))
                            with obj2.open('r+') as outfile:
                                if filetp in ('png', 'PNG'):
                                    plt.savefig(outfile, bbox_inches='tight', dpi=dpipng)
                                else:
                                    plt.savefig(outfile, bbox_inches='tight')

if __name__ == '__main__':
    cli()
