import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as mpl
import h5py

fontsize = 10

mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 8
mpl.rcParams['font.size'] = 8  # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['xtick.major.size'] = 3
mpl.rcParams['ytick.major.size'] = 3
mpl.rcParams['xtick.major.width'] = .7
mpl.rcParams['ytick.major.width'] = .7
mpl.rcParams['xtick.direction'] = 'out'
mpl.rcParams['ytick.direction'] = 'out'
mpl.rcParams['figure.titlesize'] = 8
mpl.rc('text', usetex=True)

mpl.rcParams['text.latex.preamble'] = [
    #    r'\renewcommand{\familydefault}{\sfdefault}',
    #    r'\usepackage[scaled=1]{helvet}',
    r'\usepackage[helvet]{sfmath}',
    #    r'\everymath={\sf}'
]



def getTrafoMatrix(kPoint, t1, t2, t3):

    orbitalMatrix = np.zeros((2, 2), dtype=float)
    orbitalMatrix[0, 0] = 2. * t2 * np.cos(kPoint)
    orbitalMatrix[0, 1] = 2. * t1 * np.cos(.5 * kPoint)
    orbitalMatrix[1, 0] = 2. * t1 * np.cos(.5 * kPoint)
    orbitalMatrix[1, 1] = 2. * t3 * np.cos(kPoint)

    eVals, eVecs = np.linalg.eigh(orbitalMatrix)

    return eVecs

def bandKPoint(kPoint, t1, t2, t3):

    orbitalMatrix = np.zeros((2, 2), dtype=float)
    orbitalMatrix[0, 0] = 2. * t2 * np.cos(kPoint)
    orbitalMatrix[0, 1] = 2. * t1 * np.cos(.5 * kPoint)
    orbitalMatrix[1, 0] = 2. * t1 * np.cos(.5 * kPoint)
    orbitalMatrix[1, 1] = 2. * t3 * np.cos(kPoint)

    U = getTrafoMatrix(kPoint, t1, t2, t3)

    bandMatrix = np.dot(np.conj(np.transpose(U)), np.dot(orbitalMatrix, U))

    return bandMatrix

def lmcKPoint(kPoint, t1, t2, t3):

    orbitalMatrix = np.zeros((2, 2), dtype=float)
    orbitalMatrix[0, 0] = -2. * t2 * np.sin(kPoint)
    orbitalMatrix[0, 1] = -t1 * np.sin(.5 * kPoint)
    orbitalMatrix[1, 0] = -t1 * np.sin(.5 * kPoint)
    orbitalMatrix[1, 1] = -2. * t3 * np.sin(kPoint)

    U = getTrafoMatrix(kPoint, t1, t2, t3)

    lmcMatrix = np.dot(np.conj(np.transpose(U)), np.dot(orbitalMatrix, U))

    return lmcMatrix



def bandKSpace(kVec, t1, t2, t3):
    bandStructure = np.zeros((len(kVec), 2), dtype=float)

    for kInd, kPoint in enumerate(kVec):
        bandStructure[kInd, 0] = bandKPoint(kPoint, t1, t2, t3)[0, 0]
        bandStructure[kInd, 1] = bandKPoint(kPoint, t1, t2, t3)[1, 1]


    return bandStructure

def lmcKSpace(kVec, t1, t2, t3):
    lmc = np.zeros((len(kVec), 2, 2, 1), dtype=float)

    for kInd, kPoint in enumerate(kVec):
        lmc[kInd, :, :, 0] = lmcKPoint(kPoint, t1, t2, t3)

        if(np.abs(lmc[kInd - 1, 0, 1, 0] - lmc[kInd, 0, 1, 0]) > 0.1):
            lmc[kInd, 0, 1, 0] = - lmc[kInd, 0, 1, 0]
        if(np.abs(lmc[kInd - 1, 1, 0, 0] - lmc[kInd, 1, 0, 0]) > 0.1):
            lmc[kInd, 1, 0, 0] = - lmc[kInd, 1, 0, 0]

    return lmc



def plotSawtoothBands(kVec, bandstructure):
    fig = plt.figure()
    fig.set_size_inches(3., 2.)
    ax = fig.add_subplot(111)



    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    ax.plot(kVec, bandstructure[:, 0], color = cmapBone(.5), lw=1.5)
    ax.plot(kVec, bandstructure[:, 1], color = cmapPink(.5), lw=1.5)

    ax.axhline(np.sqrt(2.) - 3., color = 'gray', lw = 1.)


    ax.set_xlim(0., 2. * np.pi)

    ax.set_xticks([0., np.pi / 2., np.pi, 3. * np.pi / 2])
    ax.set_xticklabels(['$0$', r'$\frac{\pi}{2}$', '$\pi$', r'$\frac{3\pi}{2}$'], fontsize = fontsize)

    #ax.set_ylim(-1.9, 1.7)

    ax.set_yticks([-1.5, -1., -0.5, 0., 0.5, 1., 1.5])
    ax.set_yticklabels(['$-1.5$', '$-1$', '$-0.5$', '$0$', '$0.5$', '$1$', '$1.5$'], fontsize = fontsize)

    ax.set_xlabel(r'$\rm{k}$')
    ax.set_ylabel(r'$E(\rm{k})$')

    ax.text(2.2, np.amin(bandstructure[:, 1]) + 0.6, "$W_1$ = {0:.2f}t".format(np.amax(bandstructure[:, 1]) - np.amin(bandstructure[:, 1])), fontsize = 10)
    ax.text(2.2, -1.2, "$W_0$ = {0:.2f}t".format(np.amax(bandstructure[:, 0]) - np.amin(bandstructure[:, 0])), fontsize = 10)


    #plt.tight_layout()
    #plt.show()
    plt.savefig('./savedPlots/sawtoothBands.png', format='png', bbox_inches='tight', dpi=600)

def plotSawtoothLMC(kVec, lmc):
    fig = plt.figure()
    # fig.set_size_inches(0.6 * 4., 0.6 * 2.2)
    ax = fig.add_subplot(111)
    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')
    ax.plot(kVec, lmc[:, 0, 0], color=cmapBone(.3))
    ax.plot(kVec, lmc[:, 1, 1], color=cmapBone(.7))
    ax.plot(kVec, lmc[:, 0, 1], color=cmapPink(.3))
    ax.plot(kVec, lmc[:, 1, 0], color=cmapPink(.7))
    plt.savefig('./savedPlots/sawtoothlmc.png', format='png', bbox_inches='tight', dpi=600)


def writeBandstructure(bandStructure, N, d):

    print("Bandstructure.shape() = {}".format(bandStructure.shape))

    if(d>=0):
        f = h5py.File("./input/sawtoothBandN{:d}d{:d}.h5".format(int(N), int(100 * d)), 'w')
        print("writing to file sawtoothBandN{:d}d{:d}.h5".format(int(N), int(100 * d)))
    else:
        f = h5py.File("./input/sawtoothBandN{:d}dM{:d}.h5".format(int(N), int(100 * np.abs(d))), 'w')
        print("writing to file sawtoothBandN{:d}dM{:d}.h5".format(int(N), int(100 * d)))


    dsetReal = f.create_dataset("Real", data=np.real(bandStructure))
    dsetImag = f.create_dataset("Imag", data=np.imag(bandStructure))
    f.close()

def writeLMC(lmc, N, d):

    print("LMC.shape() = {}".format(lmc.shape))
    filename = "./input/sawtoothLMCN{:d}d{:d}.h5".format(int(N), int(100 * d))
    f = h5py.File(filename, 'w')
    print(filename)

    dsetReal = f.create_dataset("Real", data=np.real(lmc))
    dsetImag = f.create_dataset("Imag", data=np.imag(lmc))
    f.close()


def main():
    print("Creating a beautiful sawtooth chain!")

    N = 1000
    d = 0.0

    kVec = np.linspace(0., 2. * np.pi, N, endpoint=False)
    t1 = -1.
    t2 = - 1. / np.sqrt(2.) + d
    t3 = 0.


    bandStructure = bandKSpace(kVec, t1, t2, t3)
    plotSawtoothBands(kVec, bandStructure)

    kVec = np.linspace(0., 2. * np.pi, 2 * N, endpoint=False)
    lmc = lmcKSpace(kVec, t1, t2, t3)

    for kInd in range(2 * N):
        lmc[kInd, 0, 0] = 0.
        #lmc[kInd, 0, 1] = 0
        #lmc[kInd, 1, 0] = 0
        lmc[kInd, 1, 1] = 0.

    plotSawtoothLMC(kVec, lmc)

    writeBandstructure(bandStructure, N, d)
    print("lower Band middle = {}".format(np.mean(bandStructure[:, 0])))
    print("E2 - Om = {}".format(np.mean(bandStructure[:, 1]) - 3.))
    writeLMC(lmc, N, d)

main()