import h5py
import numpy as np

import plotBandstructure as plotBands
import plotGFObject as plotG
import gapsOverBeta
import plotSigmaAsOfT
import dataLessPlots
import criticalTsAndGaps
import tcAsOfCoupling as tcC
import plotTcAsOfOmega
import plotSigmaIn2D
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors
import h5py
from mpl_toolkits.mplot3d import Axes3D
import pdwPlots
import plotKResolvedGap
import plotTcAsOfMu
import beautifulPlotMarios
import scipy.constants as const

def getDataFromFilename(filename):
    file = h5py.File(filename, 'r')
    gData = file['Real'][()] + 1j * file['Imag'][()]
    print('read in data.shape = ' + str(gData.shape))
    return gData


def resFactor(w1, Om, E2):
    return (E2 ** 2 - Om ** 2 - w1 ** 1) / \
           (((E2 + Om) ** 2 + w1 ** 2) * ((E2 - Om) ** 2 + w1 ** 2))


def plotResFactor():

    E2 = 5.
    N = 100

    w1Arr = np.linspace(0., 10., N)
    OmArr = np.linspace(0., 10., N)

    resFac = np.zeros((N, N))
    for w1Ind, w1 in enumerate(w1Arr):
        for omInd, Om in enumerate(OmArr):
            resFac[w1Ind, omInd] = resFactor(w1, Om, E2)

    X, Y = np.meshgrid(w1Arr, OmArr)


    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')

    ax1.set_xlabel("$\omega$")
    ax1.set_ylabel("$\Omega$")

    surf1Real = ax1.plot_surface(X, Y, resFac, ccount=256, rcount=256)

    plt.savefig('savedPlots/resFac.png', format='png', bbox_inches='tight', dpi=600)


def main():
    print("Hello from the plot program!")

    filenames = [
        "./../data/results/tcAsOfOmOneGapBos0.hdf5",
        "./../data/results/tcAsOfOmOneGapBos100.hdf5",
    ]

    beautifulPlotMarios.plotTcAsOfOmega(filenames)

    #exit()

    filename = "./../data/results/tcAsOfNphOneGapBos297.hdf5"
    beautifulPlotMarios.plotPhaseDiagram(filename)


    #convert betaC from 1/eV to K

    betaC = 5179.84
    betaC =  763.931
    TK = 1. / betaC * 1. / const.Boltzmann * const.e
    print("T[K] = {}".format(TK))

    exit()

    #throw together constants for Marios estimate
    vf = 1e6
    c1 = 1e-2
    Acell = 1e-20
    alpha = const.fine_structure
    hbarIneV = const.hbar / const.e
    print("e = {}".format(const.e))

    gSqr = hbarIneV**2 * vf**2 / Acell * c1 * alpha * 4. * np.pi * const.c / (4 * 1e12) / 1e-8
    print("gSqr = {}".format(gSqr))




main()
