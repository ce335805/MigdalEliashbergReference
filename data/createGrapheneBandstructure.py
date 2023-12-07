import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
import h5py
import scipy.optimize
from matplotlib.colors import ListedColormap
import scipy.constants as consts

# real space lattice vectors
a1 = np.array([np.sqrt(3.) / 2., 3. / 2.])
a2 = np.array([np.sqrt(3.), 0.])

# reciprocal lattice vectors
b1 = 4. * np.pi / (3. * np.sqrt(3.)) * np.array([0., np.sqrt(3.)])
b2 = 4. * np.pi / (3. * np.sqrt(3.)) * np.array([3. / 2., -np.sqrt(3.) / 2.])

N = 50
writeFrac = 1. / 500.
delta =0.014
tHop = 2.7
tDash = -0. * tHop
mu = 0.
dist = 250. / 0.142

# nearest neighbour vectors
nn1 = np.array([np.sqrt(3.) / 2., 1. / 2.])
nn2 = np.array([-np.sqrt(3.) / 2., 1. / 2.])
nn3 = np.array([0., -1.])

def orbitalMatrix(kPoint):

    diagElem = tDash * (np.exp(1j * np.dot(kPoint, a1)) + np.exp(1j * np.dot(kPoint, a2)) + np.exp(1j * np.dot(kPoint, a1 + a2)))
    diagElem = diagElem + np.conj(diagElem)

    orbMat = np.array([
        [delta / 2. + diagElem,
         tHop * (np.exp(1j * np.dot(kPoint, nn1)) + np.exp(1j * np.dot(kPoint, nn2)) + np.exp(
             1j * np.dot(kPoint, nn3)))],
        [tHop * (np.exp(-1j * np.dot(kPoint, nn1)) + np.exp(-1j * np.dot(kPoint, nn2)) + np.exp(
            -1j * np.dot(kPoint, nn3))),
         -delta / 2. + np.conj(diagElem)]
    ])

    return orbMat

def bandStructure(kPoint):
    orbMat = orbitalMatrix(kPoint)

    val, vec = np.linalg.eigh(orbMat)

    if(val[0] > val[1]):
        val = np.flip(val)

    return val

def dkxOrbitalMatrix(kPoint):
    
    diagElem = tDash * (1j * a1[0] * np.exp(1j * np.dot(kPoint, a1)) + 
                        1j * a2[0] * np.exp(1j * np.dot(kPoint, a2)) + 
                        1j * (a1[0] + a2[0]) * np.exp(1j * np.dot(kPoint, a1 + a2)))
    diagElem = diagElem + np.conj(diagElem)
    
    orbMat = np.array([
        [diagElem,
         tHop * ( 1j * nn1[0] * np.exp(1j * np.dot(kPoint, nn1)) + 1j * nn2[0] * np.exp(1j * np.dot(kPoint, nn2)) + 1j * nn3[0] * np.exp(
             1j * np.dot(kPoint, nn3)))],
        [tHop * ( -1j * nn1[0] * np.exp(-1j * np.dot(kPoint, nn1)) - 1j * nn2[0] * np.exp(-1j * np.dot(kPoint, nn2)) - 1j * nn3[0] * np.exp(
            -1j * np.dot(kPoint, nn3))),
         np.conj(diagElem)]
    ])

    return orbMat

def dkyOrbitalMatrix(kPoint):
    diagElem = tDash * (1j * a1[1] * np.exp(1j * np.dot(kPoint, a1)) +
                        1j * a2[1] * np.exp(1j * np.dot(kPoint, a2)) +
                        1j * (a1[1] + a2[1]) * np.exp(1j * np.dot(kPoint, a1 + a2)))
    diagElem = diagElem + np.conj(diagElem)

    orbMat = np.array([
        [diagElem,
         tHop * ( 1j * nn1[1] * np.exp(1j * np.dot(kPoint, nn1)) + 1j * nn2[1] * np.exp(1j * np.dot(kPoint, nn2)) + 1j * nn3[1] * np.exp(
             1j * np.dot(kPoint, nn3)))],
        [tHop * ( -1j * nn1[1] * np.exp(-1j * np.dot(kPoint, nn1)) - 1j * nn2[1] * np.exp(-1j * np.dot(kPoint, nn2)) - 1j * nn3[1] * np.exp(
            -1j * np.dot(kPoint, nn3))),
         diagElem]
    ])

    return orbMat

def lmc(kPoint):
    orbMat = orbitalMatrix(kPoint)
    val, vec = np.linalg.eigh(orbMat)

    lmcXMat = dkxOrbitalMatrix(kPoint)
    lmcYMat = dkyOrbitalMatrix(kPoint)

    lmcX = np.matmul(np.transpose(np.conj(vec)), np.matmul(lmcXMat, vec))
    lmcY = np.matmul(np.transpose(np.conj(vec)), np.matmul(lmcYMat, vec))

    return (lmcX, lmcY)

def calcCouplingSinglePoint(kPoint):
    coup = np.sqrt(np.dot(kPoint, kPoint)) * np.exp(- dist * np.sqrt(np.dot(kPoint, kPoint)))
    return coup

def calcBoxCouplingAroundK(kPoint, fraction, aOverD):
    KPoint = b1 / 3. + 2. * b2 / 3.
    abs = np.dot(kPoint - KPoint, kPoint - KPoint)
    if(np.sqrt(abs) < aOverD):
        return 1.
    return 0.


def calcBoxFunc(kPoint):
    if (np.sqrt(np.dot(kPoint, kPoint)) > 0.003):
        return 0
    else:
        return 1

def plotBandstructure(vec1, vec2, bandstructure, boxCoupling):

    coordMat = np.zeros((N, N, 2))

    for ind1, v1 in enumerate(vec1):
        for ind2, v2 in enumerate(vec2):
            coordMat[ind1, ind2, :] = v1 + v2

    coordMatX = coordMat[:, :, 0]
    coordMatY = coordMat[:, :, 1]

    print("coordMatX.shape = {}".format(coordMatX.shape))
    print("coordMatY.shape = {}".format(coordMatY.shape))
    print("bandstructure.shape = {}".format(bandstructure.shape))

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    fig.set_size_inches(4., 3.)
    plt.tight_layout()

    bandstructure = np.reshape(bandstructure, (N, N, 2))
    boxCoupling = np.reshape(boxCoupling, (N, N))

    cmapCoolWarm = cm.coolwarm
    cmapCoolWarmTrans = cmapCoolWarm(np.arange(cmapCoolWarm.N))
    cmapCoolWarmTrans[:, -1] = 0.5
    cmapCoolWarmTrans = ListedColormap(cmapCoolWarmTrans)

    #surfBox = ax1.plot_surface(coordMatX, coordMatY, np.real(boxCoupling[:, :]), cmap=cm.viridis, ccount=256, rcount=256)

    surf1Real = ax1.plot_surface(coordMatX, coordMatY, np.real(bandstructure[:, :, 0]), cmap=cm.coolwarm, ccount=256, rcount=256)
    surf2Real = ax1.plot_surface(coordMatX, coordMatY, np.real(bandstructure[:, :, 1]), cmap=cm.coolwarm_r, ccount=256, rcount=256)

    ax1.set_xlabel("$k_x$", fontsize = 10)
    ax1.set_ylabel("$k_y$", fontsize = 10)
    ax1.set_zlabel(r"$\varepsilon(k)$", fontsize = 10)

    ax1.view_init(elev=10, azim=30)
    ax1.dist = 11
    #plt.show()
    plt.savefig('savedPlots/GrapheneBands.png', format='png', bbox_inches='tight', dpi = 600)

def plotBandstructureFancy(vec1, vec2, bandstructure, boxCoupling):

### for plot in paper used d=25nm

    coordMat = np.zeros((N, N, 2))

    for ind1, v1 in enumerate(vec1):
        for ind2, v2 in enumerate(vec2):
            coordMat[ind1, ind2, :] = v1 + v2

    coordMatX = coordMat[:, :, 0]
    coordMatY = coordMat[:, :, 1]

    print("coordMatX.shape = {}".format(coordMatX.shape))
    print("coordMatY.shape = {}".format(coordMatY.shape))
    print("bandstructure.shape = {}".format(bandstructure.shape))

    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    fig.set_size_inches(4., 3.)
    plt.tight_layout()

    bandstructure = np.reshape(bandstructure, (N, N, 2))
    boxCoupling = np.reshape(boxCoupling, (N, N))

    cmapCoolWarm = cm.coolwarm
    cmapCoolWarmTrans = cmapCoolWarm(np.arange(cmapCoolWarm.N))
    cmapCoolWarmTrans[:, -1] = 0.4
    cmapCoolWarmTrans[110:, -1] = 0.
    cmapCoolWarmTrans = ListedColormap(cmapCoolWarmTrans)

    cmapCoolWarm_r = cm.coolwarm_r
    cmapCoolWarmTrans_r = cmapCoolWarm_r(np.arange(cmapCoolWarm_r.N))
    cmapCoolWarmTrans_r[:, -1] = 0.4
    cmapCoolWarmTrans_r[: 256 - 110, -1] = 0.
    cmapCoolWarmTrans_r = ListedColormap(cmapCoolWarmTrans_r)

    viridis = cm.get_cmap('viridis', 256)
    newcolors = viridis(np.linspace(0, 1, 256))
    opaque = np.array([1., 1., 1., 0.])
    newcolors[:150, :] = opaque
    #newcolors[:-25, :] = opaque
    newcmp = ListedColormap(newcolors)

    #surfBox = ax1.plot_surface(coordMatX, coordMatY, np.real(boxCoupling[:, :]), cmap=cm.viridis, ccount=256, rcount=256)

    surf1Real = ax1.plot_surface(coordMatX, coordMatY, np.real(bandstructure[:, :, 0]), cmap=cmapCoolWarmTrans_r, ccount=256, rcount=256)
    surf2Real = ax1.plot_surface(coordMatX, coordMatY, np.real(bandstructure[:, :, 1]), cmap=cmapCoolWarmTrans, ccount=256, rcount=256)

    ax1.axhline(-0.01, 0.45, 0.6, color = "red", lw = 0.8, linestyle = '--', zorder = 666)
    ax1.axhline(+0.005, 0.5, 0.6, color = "red", lw = 0.8, linestyle = '-', zorder = 666)

    ax1.text(0., 1.9, 7.7, r"$\mu$", transform = ax1.transAxes, fontsize = 12)
    ax1.text(0., 2.6, 8.1, r"$\Omega_{\rm s}$", transform = ax1.transAxes, color = '#666666', fontsize = 12)

    ax1.arrow(0.55, 0.46, 0., 0.05, head_length = 0.01, head_width = 0.01, transform = ax1.transAxes, color = '#666666')

    ax1.set_xticks([])
    ax1.set_yticks([])
    ax1.set_zticks([])
    ax1.set_zlim(-0.07, 0.07)

    #ax1.set_xlabel("$k_x$", fontsize = 10)
    #ax1.set_ylabel("$k_y$", fontsize = 10)
    #ax1.set_zlabel(r"$\varepsilon(k)$", fontsize = 10)

    #ax1.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))

    ax1.xaxis.pane.fill = False
    ax1.yaxis.pane.fill = False
    ax1.zaxis.pane.fill = False
    ax1.xaxis.pane.set_edgecolor('w')
    ax1.yaxis.pane.set_edgecolor('w')
    ax1.zaxis.pane.set_edgecolor('w')
    ax1.grid(False)

    ax1.set_axis_off()

    for axis in ['top', 'right', 'bottom', 'left']:
        ax1.spines[axis].set_linewidth(0.0)

    ax1.view_init(elev=0, azim=0)
    ax1.dist = 11
    #plt.show()
    plt.savefig('savedPlots/GrapheneBandsFancy.png', format='png', bbox_inches='tight', dpi = 600)



def createKPath():
    pathN = int(N / writeFrac)
    path = np.zeros((pathN + 1, 2))
    choppedB1 = np.zeros((pathN, 2))
    choppedB2 = np.zeros((pathN, 2))
    for ind in range(pathN):
        choppedB1[ind, :] = b1 / pathN * ind
        choppedB2[ind, :] = b2 / pathN * ind

    for ind in range(pathN // 3 + 1):
        path[ind, :] = 2. * b2 / pathN * ind + b1 / pathN * ind

    for ind in np.arange(pathN // 6):
        path[ind + pathN // 3 + 1, :] = path[pathN // 3] - 2. * b1 / pathN * (ind + 1) - b2 / pathN * (ind + 1)

    for ind in range(pathN // 2):
        path[ind + pathN // 3 + pathN // 6 + 1, :] = path[pathN // 3 + pathN // 6] - b2 / pathN * (ind + 1)

    return path



def plotBandstructureOnPath():
    pathN = int(N / writeFrac)
    kPath = createKPath()
    bandsOnPath = np.zeros((pathN + 1, 2))
    lmcXOnPath = np.zeros((pathN + 1, 2, 2), dtype=complex)
    lmcYOnPath = np.zeros((pathN + 1, 2, 2), dtype=complex)
    couplingsRelK = np.zeros(pathN + 1)
    boxRelK = np.zeros(pathN + 1)
    for kInd, kPoint in enumerate(kPath):
        bandsOnPath[kInd, :] = bandStructure(kPoint)
        lmcXOnPath[kInd, :] = lmc(kPoint)[0]
        lmcYOnPath[kInd, :] = lmc(kPoint)[1]
        couplingsRelK[kInd] = calcCouplingSinglePoint(kPoint - kPath[pathN//3])
        boxRelK[kInd] = calcBoxFunc(kPoint - kPath[pathN//3])

    lmcOnPath = lmcXOnPath + lmcYOnPath

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(4., 3.)
    plt.tight_layout()
    ax2 = ax.twinx()

    ax.plot(np.arange(pathN + 1), bandsOnPath[:, 0], color = 'teal', linewidth = 1.)
    ax.plot(np.arange(pathN + 1), bandsOnPath[:, 1], color = 'goldenrod', linewidth = 1.)

    ax.axhline(mu, color = 'red', linewidth = 0.3, linestyle = '-')

    #ax2.plot(np.arange(N + 1), lmcOnPath[:, 0, 0], color = 'teal', linewidth = 1.)
    #ax2.plot(np.arange(N + 1), lmcOnPath[:, 1, 1], color = 'goldenrod', linewidth = 1.)
    ax2.plot(np.arange(pathN + 1), np.abs(lmcOnPath[:, 0, 1]), color = 'gray', linewidth = 1.)
    ax2.plot(np.arange(pathN + 1), np.abs(lmcOnPath[:, 1, 0]), color = 'gray', linewidth = 1.)
    ax2.plot(np.arange(pathN + 1), 50. * np.abs(couplingsRelK), color = 'mediumseagreen', linewidth = 1.)
    ax2.plot(np.arange(pathN + 1), boxRelK, color = 'goldenrod', linewidth = 1.)

    ax.axvline(pathN // 3, linestyle='--', linewidth=0.5, color='gray')
    ax.axvline(pathN // 2, linestyle='--', linewidth=0.5, color='gray')
    #ax.axhline(0., linewidth=0.5, color='black')
    ax.set_xticks([0, pathN // 3, pathN // 2, pathN])
    ax.set_xticklabels([r"$\Gamma$", r"$K$", r"$M$", r"$\Gamma$"])

    ax.set_xlim(pathN//3 - pathN * writeFrac * 2., pathN//3 + pathN * writeFrac * 2.)
    ax.set_ylim(-.1, .1)
    ax2.set_ylim(0.)

    plt.savefig('savedPlots/grapheneBandsOnPath.png', format='png', bbox_inches='tight', dpi=600)

def effectiveDensityOfStates():
    bands = np.zeros((N * N, 2), dtype='complex')
    coupFromK = np.zeros((N * N), dtype='complex')

    choppedB1 = np.zeros((N, 2))
    choppedB2 = np.zeros((N, 2))
    for ind in range(N):
        choppedB1[ind, :] = b1 / N * ind
        choppedB2[ind, :] = b2 / N * ind

    for ind1, cB1 in enumerate(choppedB1):
        for ind2, cB2 in enumerate(choppedB2):
            bands[ind1 * N + ind2, :] = bandStructure(cB1 + cB2)
            #coupFromK[ind1 * N + ind2] = calcCouplingSinglePoint(cB1 + cB2 - b1 / 3. - 2. * b2 / 3.)
            coupFromK[ind1 * N + ind2] = 1.



    EArr = np.linspace(-10., 10., 200, endpoint=True)

    DKE = np.zeros(len(EArr))

    eta = 5. * 1e-2
    for eInd, eVal in enumerate(EArr):
        deltaFunc = 1. / (2. * np.pi) * eta / ((eVal - bands) ** 2 + (eta / 2.) ** 2)
        DKE[eInd] += 1 / N ** 2 * np.sum(np.real(coupFromK * deltaFunc[:, 0]))
        DKE[eInd] += 1 / N ** 2 * np.sum(np.real(coupFromK * deltaFunc[:, 1]))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(4., 3.)
    plt.tight_layout()

    ax.plot(EArr, DKE)

    ax.axvline(mu, linewidth=0.5, color='black')
    ax.plot(np.linspace(1.2 - 2.7, 1.2, 100), np.ones(100) * 0.5, linewidth = 0.5, color = 'red')
    #ax.set_xlim(-1, 1)
    ax.set_ylim(0., 4.)

    plt.savefig('savedPlots/Deff.png', format='png', bbox_inches='tight', dpi=600)



def writeBandstructureToFile(fraction, aOverD, d):

    bands = np.zeros((N * N, 2), dtype='complex')
    coupling = np.zeros((N * N), dtype='complex')

    choppedB1 = np.zeros((N, 2))
    choppedB2 = np.zeros((N, 2))
    for ind in range(N):
        #choppedB1[ind, :] = b1 / N * ind
        #choppedB2[ind, :] = b2 / N * ind
        #choppedB1[ind, :] = b1 / N * writeFrac * ind + b1 / 3. - b1 * writeFrac / 2.
        #choppedB2[ind, :] = b2 / N * writeFrac * ind + 2. * b2 / 3. - b2 * writeFrac / 2.
        choppedB1[ind, :] = b1 / N * 1. / fraction * ind + b1 / 3. - b1 / fraction / 2.
        choppedB2[ind, :] = b2 / N * 1. / fraction * ind + 2. * b2 / 3. - b2 / fraction / 2.

    print(b1 / 4. / fraction + b1 / 3. + b2 / 4. / fraction + 2. * b2 / 3.)

    for ind1, cB1 in enumerate(choppedB1):
        for ind2, cB2 in enumerate(choppedB2):
            bands[ind1 * N + ind2, :] = bandStructure(cB1 + cB2)
            coupling[ind1 * N + ind2] = calcBoxCouplingAroundK(cB1 + cB2, fraction, aOverD)

    print("max Band = {}".format(np.amax(bands)))
    print("min Band = {}".format(np.amin(bands)))

    bandsTimesCoupling = np.zeros(bands.shape)
    for ind in np.arange(bands.shape[0]):
        bandsTimesCoupling[ind, 0] = bands[ind, 0]# * coupling[ind]
        bandsTimesCoupling[ind, 1] = bands[ind, 1]# * coupling[ind]
        #bandsTimesCoupling[ind, 0] = coupling[ind]
        #bandsTimesCoupling[ind, 1] = coupling[ind]

    print("bandsMin[0] = {}".format(np.amin(bands[:, 0])))
    print("bandsMax[0] = {}".format(np.amax(bands[:, 0])))

    print("bandsMin[1] = {}".format(np.amin(bands[:, 1])))
    print("bandsMax[1] = {}".format(np.amax(bands[:, 1])))

    print("grapheneBands.shape() = {}".format(bands.shape))
    plotBandstructure(choppedB1, choppedB2, bandsTimesCoupling, coupling)
    #plotBandstructureFancy(choppedB1, choppedB2, bandsTimesCoupling, coupling)

    avGap = 1 / (N * N) * (np.real(np.sum(bands[:, 1]) - np.sum(bands[:, 0])))
    print("average band-gap = {}".format(avGap))

    #calculate approximate mu
    muApprox = np.array([mu])
    print("approx mu = {}".format(muApprox))
    print("approx mu.shape = {}".format(muApprox.shape))

    #calculate approximate bandwidth
    wApprox = np.amax(bands[:, 0]) - np.amin(bands[:, 0])
    print("W approx = {}".format(np.real(wApprox)))

    filename = "input/grapheneN{}D{}d{}.h5".format(int(N), int(delta * 1000), int(d))
    #filename = "input/flatN{}D{}.h5".format(int(N), int(delta * 10))
    f = h5py.File(filename, 'w')
    print("Writing to file: " + filename)
    dsetReal = f.create_dataset("Real", data=np.real(bands))
    dsetImag = f.create_dataset("Imag", data=np.imag(bands))

    dsetMu = f.create_dataset("MuApprox", data=muApprox)
    f.close()

    filename = "input/boxCouplingN{}d{}.h5".format(int(N), int(d))
    f = h5py.File(filename, 'w')
    print("Writing to file: " + filename)
    dsetReal = f.create_dataset("Real", data=np.real(coupling))
    dsetImag = f.create_dataset("Imag", data=np.imag(coupling))
    f.close()


def estimateHeating(d, nBos):

    a = 0.142 * 1e-9
    aSTO = 0.3905 * 1e-9
    xi = d
    #cP = 5 * 1e-3
    T = 15
    cP = 0.75 * T**3 * 1e-4
    energyPerBoson = 3.9 * 4.136 * 1e-3 * consts.e

    AGraph = np.abs(np.dot(a1, a2) * a**2)
    fracBZ = np.abs(np.abs(1. / (np.dot(b1, b2) / a**2) / d**2 * np.pi))

    print("fraction of BZ couples: {}".format(fracBZ))
    print("relative volumes: {}".format(AGraph * xi / aSTO**3))
    print("fraction and volume = {}".format(fracBZ *AGraph * xi / aSTO**3))

    dT = energyPerBoson * nBos * fracBZ * aSTO**3 / AGraph / xi * consts.Avogadro / cP
    #dT = energyPerBoson * nBos * fracBZ / aSTO**3 * AGraph * xi * consts.Avogadro / cP

    print("expected change in T: {}".format(dT))

    ### estimate heating assuming all heat goes into graphene

    cPGraph = 0.24#Heat capacity at 10K in J / mol / K
    dTGraph = energyPerBoson * nBos * fracBZ * consts.Avogadro / cPGraph
    print("Temperature change graphene = {}".format(dTGraph))

def main():
    print("Producing beautiful gaphene bands!")


    print("Estimate energy connected by q: {}".format(consts.hbar * 1e14 / consts.e))

    print("fraction of BZ = {}".format(np.floor(dist * 4. * np.pi / 3.)))

    d = 5.#nm
    #initial 0.5 just to be sure I have enough points
    fractionFromD = 0.5 * 1. / 2. * d / 0.142 * 4. * np.pi / 3.

    estimateHeating(d * 1e-9, 1.)

    #exit()

    print("BZ fraction considered: {}".format(fractionFromD))

    writeBandstructureToFile(fractionFromD, 0.142 / d, d)

    #effectiveDensityOfStates()

    vf = 1e6
    a = 1e-9 * 0.142
    wT = 1.26 * 1e12 * 2. * np.pi
    wL = 5.1 * 1e12 * 2. * np.pi
    epsInf = 6.3
    epsHBN = 5
    hbarIneV = consts.hbar / consts.e


    wLim = np.sqrt((epsHBN * wT**2 + epsInf * wL**2) / (epsInf + epsHBN))
    print("wLim = {}THz".format(wLim * 1e-12))
    print("wLim in eV = {}".format(wLim * hbarIneV))
    fw = ( 4. * wLim**2 * (epsInf + epsHBN)) / (wLim**2 - wT**2)
    fw = 1. / np.sqrt(fw)
    print("fw = {}".format(fw))

    fac1 = np.sqrt(np.pi * consts.fine_structure) * np.sqrt(2.) / np.sqrt(3 * np.sqrt(3))
    fac2 = hbarIneV * vf * np.sqrt(consts.c) / (a * np.sqrt(wLim))
    gTimesSqrtD = fac1 * fac2 * fw
    print("My g * Sqrt(d) = {}".format(gTimesSqrtD))
    print("My g = {}eV".format(gTimesSqrtD * np.sqrt(1. / (d * 1e-9))))


    print("g^2 = {}eV^2".format((gTimesSqrtD * np.sqrt(1. / (d * 1e-9)))**2))
    print("g^2 = {}eV^2 nm^2".format((gTimesSqrtD * np.sqrt(1. / (d * 1e-9)))**2 * 3. * np.sqrt(3.) / 2. * a**2 / (2. * np.pi)))

    estimate2 = 1. / (2. * np.pi) * np.pi * consts.fine_structure * fw ** 2 * hbarIneV ** 2 * vf ** 2 * consts.c / (d * 1e-9 * wLim)

    print("estimate2 = {}".format(estimate2 * 1e18))
main()
