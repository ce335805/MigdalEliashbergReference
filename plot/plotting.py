import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
import matplotlib.patches as patches
import matplotlib as mpl
import matplotlib.cm as cm
import matplotlib.colors
import h5py
from matplotlib import gridspec
from matplotlib.patches import ConnectionPatch
import matplotlib.image as mpimg
from matplotlib import pyplot
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon

fontsize = 8

mpl.rcParams['font.family'] = 'Helvetica'
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['lines.markersize'] = 8
mpl.rcParams['font.size'] = 8  # <-- change fonsize globally
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['axes.titlesize'] = 8
mpl.rcParams['axes.labelsize'] = 8
mpl.rcParams['xtick.major.size'] = 3
mpl.rcParams['ytick.major.size'] = 3
mpl.rcParams['xtick.major.width'] = .5
mpl.rcParams['ytick.major.width'] = .5
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


def getCandBFromFile(fileName):
    file = h5py.File(fileName, 'r')
    couplings = file['Couplings'][()]
    betaCs = file['BetaCs'][()]

    return (couplings, betaCs)

def computeDensityOfStatesFromBands(bands, N):
    EArr = np.linspace(-1.7, 1., 200, endpoint=True)

    rho = np.zeros(len(EArr))

    eta = 5. * 1e-3
    for eInd, eVal in enumerate(EArr):
        deltaFunc = 1. / (2. * np.pi) * eta / ((eVal - bands) ** 2 + (eta / 2.) ** 2)
        rho[eInd] += 1 / N ** 2 * np.sum(np.real(deltaFunc[:, 0]))
        rho[eInd] += 1 / N ** 2 * np.sum(np.real(deltaFunc[:, 1]))
    return rho, EArr

def plotTcAsOfOmega(fileNames):

    nBosArr = np.array([0., 0.1])

    omegasLen, _ = getCandBFromFile(fileNames[0])
    nOmega = len(omegasLen)

    omegaArrs = np.zeros((len(fileNames), nOmega))
    betaCArrs = np.zeros((len(fileNames), nOmega))

    for fileInd, fileName in enumerate(fileNames):
        omegas, betasCs = getCandBFromFile(fileName)
        omegaArrs[fileInd, :] = omegas
        betaCArrs[fileInd, :] = betasCs

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(2.5, 1.8)

    cmapBone = cm.get_cmap('bone')
    cmapPink = cm.get_cmap('pink')

    ax.plot(omegaArrs[0, :], 1. / betaCArrs[0, :], color=cmapPink(0.3), linestyle='-', linewidth=1.5,
                label=r"$n_{\rm Drive}{=}0$", zorder = 6)

    ax.plot(omegaArrs[1, :], 1. / betaCArrs[1, :], color=cmapPink(0.5), linestyle='-', linewidth=1.5,
                label=r"$n_{\rm Drive}{=}1$")

    ax.plot(omegaArrs[2, :], 1. / betaCArrs[2, :], color=cmapBone(0.5), linestyle='--', dashes=[4, 4], linewidth=1.,
                label=r"$n_{\rm Drive} {=} 1 , \, \rm{Eq.}12$")

    ax.set_xlabel(r"$\Omega[\xi_d]$", fontsize=8)
    ax.set_ylabel(r"$T_c[\xi_d]$", fontsize=8, labelpad=-3)

    #ax.set_yscale('log')

    ax.set_ylim(- 3. * 1e-5, 0.003)
    ax.set_xlim(2.5, 3.07)

    ax.set_xticks([2.7, 3.])
    ax.set_xticklabels([r'$0.9$', r'$1$'], fontsize=8)
    ax.set_yticks([0., 0.001 * 3])
    ax.set_yticklabels([r'$0$', r'$0.001$'], fontsize=8)

    #
    #boxProps = dict(boxstyle='square', facecolor='white', alpha=0., linewidth=0., fill=True, pad=0.1)
    #ax.text(0.61, 1.02, r"$\Omega=\Delta$", fontsize=8, bbox=boxProps, transform = ax.transAxes)

    #ax.text(1.53, 0.004, r"$\rm{Eq. \, X}$", fontsize = 8)
    #ax.arrow(1.525, 0.0041, -0.03, -0.0001, lw = 0.1, head_width = 0., head_length=0., width=0.00002, color = 'black')

    handles, labels = ax.get_legend_handles_labels()
    legend = ax.legend(handles[:-1], labels[:-1], fontsize=fontsize - 2, loc='upper left', bbox_to_anchor=(0., 1.), edgecolor='black', ncol=1, handlelength = 2.)
    legend.get_frame().set_alpha(0.)
    legend.get_frame().set_boxstyle('Square', pad=0.1)
    legend.get_frame().set_linewidth(0.0)

    legend2 = ax.legend([handles[-1]], [labels[-1]], fontsize=fontsize - 2, loc='upper right', bbox_to_anchor=(0.87, 1.), edgecolor='black', ncol=1, handlelength = 2.)
    legend2.get_frame().set_alpha(0.)
    legend2.get_frame().set_boxstyle('Square', pad=0.1)
    legend2.get_frame().set_linewidth(0.0)

    pyplot.gca().add_artist(legend)
    ax.axvline(3., color = "gray", linestyle = '-',lw = 0.6, zorder = 666)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)

    kArr = np.linspace(-np.pi, np.pi, 1000, endpoint=True)
    eps1 = -np.sqrt(2.) * (np.cos(kArr) + 1)
    eps2 = np.ones(kArr.shape) * np.sqrt(2.)

    inax = fig.add_axes([0.25, 0.25, 0.4, 0.4], anchor='NE', zorder=666)

    inax.plot(kArr, eps1, color = 'gray', lw = 0.7)
    inax.plot(kArr, eps2, color = 'gray', lw = 0.7)

    inax.axhline(-1.5857864376269044, linestyle = '--', color = 'darkolivegreen', lw = 0.5)

    inax.set_xlabel(r"$k$", fontsize = 8)
    #inax.set_ylabel(r"$\varepsilon$", fontsize = 8)

    inax.set_xlim(-np.pi, np.pi)

    inax.set_yticks([-1.5857864376269044, np.sqrt(2.)])
    inax.set_yticklabels([r"$\mu$", r"$\varepsilon_d$"], fontsize = 8)
    inax.set_xticks([])

    inax.arrow(-1.47, -1.5857864376269044, 0., 2.9, lw = 0.6, color = "peru", head_width = 0.2, length_includes_head = True)

    inax.text(-1., -0.3, r"$\xi_d$", color = 'peru', fontsize = 8)

    for axis in ['top', 'bottom', 'left', 'right']:
        inax.spines[axis].set_linewidth(0.5)

    ax.text(-0.25, 1.07, r"$\rm{b.)}$", transform = ax.transAxes, fontsize = 10)

    # plt.tight_layout()
    # plt.show()
    plt.savefig('./savedPlots/tcAsOfOm.png', format='png', bbox_inches='tight', dpi=600)



def plotPhaseDiagram(fileName):

    nPh, betaC = getCandBFromFile(fileName)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    fig.set_size_inches(2.5, 1.8)

    ax.plot(nPh, 1. / betaC, color='black', linestyle='-', linewidth=0.5)
    ax.fill_between(nPh, 1. / betaC, color = 'teal', alpha = 0.4)
    #ax.fill_between(nPh, 1. / betaC, np.ones(betaC.shape) * 1., color = 'peru', alpha = 0.4)

    ax.set_xlabel(r"$\rm{Driving \,\, Strength}$" + r" $(N_{\rm Drive})$", fontsize=8)
    ax.set_ylabel(r"$T \, [\xi_d]$", fontsize=8)

    ax.set_yticks([0., 0.002 * 3, 0.004 * 3])
    ax.set_yticklabels([r"$0$", r"$0.002$", r"$0.004$"])

    ax.set_ylim(0., 1. / betaC[-1] + 0.005)
    ax.set_xlim(0., nPh[-1])

    ax.text(3, 0.0005, r"$\Omega = 0.99 \xi_d$", fontsize = 10)
    ax.text(0.75, 0.004, r"$\rm{NS}$", fontsize = 12)
    ax.text(3., 0.003, r"$\rm{SC}$", fontsize = 12)

    for axis in ['top', 'bottom', 'left', 'right']:
        ax.spines[axis].set_linewidth(0.5)


    ### Inset driven occupation
    xArr = np.linspace(0., 2., 1000)

    yArr = 1. / (np.exp(xArr) - 1.)
    gam = 0.2
    x_0 = 1.
    strength = 2.
    lorentz = strength * 1. / np.pi * 0.5 * gam / ((xArr - x_0)**2 + 0.25 * gam**2)

    inax1 = fig.add_axes([0.16, 0.62, 0.2, 0.2], anchor='NE', zorder= 666)
    inax2 = fig.add_axes([0.57, 0.62, 0.2, 0.2], anchor='NE', zorder= 666)

    inax1.plot(xArr, yArr, color = 'black', lw = 1.)
    inax1.set_ylim(0., 10.)
    inax1.set_xlim(0., 2.)
    inax1.set_yticks([])
    inax1.set_xticks([1])
    inax1.set_xticklabels([r"$\Omega$"], fontsize=7)

    inax1.text(0.4, 9., r"$n_{\rm B}$", fontsize = 6)
    inax1.arrow(1.7, 5., 1., 0.,fc='k', ec='k', lw = 0.6, head_width = 0.9, head_length = 0.2, overhang = 0.5,
             length_includes_head= True, clip_on = False)

    xmin, xmax = inax1.get_xlim()
    ymin, ymax = inax1.get_ylim()

    # get width and height of axes object to compute
    # matching arrowhead length and width
    dps = fig.dpi_scale_trans.inverted()
    bbox = ax.get_window_extent().transformed(dps)
    width, height = bbox.width, bbox.height

    hw = 1./20.*(ymax-ymin)
    hl = 1./20.*(xmax-xmin)
    lw = 1. # axis line width
    ohg = 0.3 # arrow overhang

    # compute matching arrowhead length and width
    yhw = hw/(ymax-ymin)*(xmax-xmin)* height/width
    yhl = hl/(xmax-xmin)*(ymax-ymin)* width/height


    inax1.arrow(0., 0., 2.2, 0., fc='k', ec='k', lw = 0.5,
             head_width=hw, head_length=hl, overhang = ohg,
             length_includes_head= True, clip_on = False)

    inax1.arrow(0., 0., 0., 11., fc='k', ec='k', lw = 0.5,
             head_width=yhw, head_length= yhl, overhang = ohg,
             length_includes_head= True, clip_on = False)


    inax2.arrow(0., 0., 2.2, 0., fc='k', ec='k', lw = 0.5,
             head_width=hw, head_length=hl, overhang = ohg,
             length_includes_head= True, clip_on = False)

    inax2.arrow(0., 0., 0., 11., fc='k', ec='k', lw = 0.5,
             head_width=yhw, head_length= yhl, overhang = ohg,
             length_includes_head= True, clip_on = False)

    gradient_fill(xArr, lorentz, "peru", inax2)

    inax2.plot(xArr, yArr + lorentz, color = 'black', lw = 1.)
    inax2.axhline(6.5, lw = 0.5, color = 'gray', linestyle = '--')
    inax2.set_ylim(0., 10.)
    inax2.set_xlim(0., 2.)
    inax2.set_yticks([6.5])
    inax2.set_yticklabels([r"$n_{\rm Drive}$"], fontsize = 6)
    inax2.set_xticks([1])
    inax2.set_xticklabels([r"$\Omega$"], fontsize = 7)


    inax2.text(0.4, 9., r"$n_{\rm B} + n_{\rm Drive}$", fontsize = 6)


    for axis in ['top', 'bottom', 'left', 'right']:
        inax1.spines[axis].set_linewidth(0.5)
        inax2.spines[axis].set_linewidth(0.5)
    for axis in ['top', 'right']:
        inax1.spines[axis].set_linewidth(0.0)
        inax2.spines[axis].set_linewidth(0.0)

    ###Add Sawtooth chain imagee
    #img = mpimg.imread('./SawtoothChain.png')
    #imax = fig.add_axes([0.13, 0.39, 0.48, 0.48], anchor='NE', zorder= 666)
    #imax.imshow(img)
#
    #for axis in ['top', 'bottom', 'left', 'right']:
    #    imax.spines[axis].set_linewidth(0.0)
    #imax.set_xticks([])
    #imax.set_yticks([])

    ax.text(-0.28, 1.07, r"$\rm{a.)}$", transform = ax.transAxes, fontsize = 10)

    plt.savefig('./savedPlots/PhaseDiagram.png', format='png', bbox_inches='tight', dpi=600)


def gradient_fill(x, y, fill_color=None, ax=None, **kwargs):
    """
    Plot a line with a linear alpha gradient filled beneath it.

    Parameters
    ----------
    x, y : array-like
        The data values of the line.
    fill_color : a matplotlib color specifier (string, tuple) or None
        The color for the fill. If None, the color of the line will be used.
    ax : a matplotlib Axes instance
        The axes to plot on. If None, the current pyplot axes will be used.
    Additional arguments are passed on to matplotlib's ``plot`` function.

    Returns
    -------
    line : a Line2D instance
        The line plotted.
    im : an AxesImage instance
        The transparent gradient clipped to just the area beneath the curve.
    """
    if ax is None:
        ax = plt.gca()

    line, = ax.plot(x, y, lw = 0.)
    if fill_color is None:
        fill_color = line.get_color()

    zorder = line.get_zorder()
    alpha = line.get_alpha()
    alpha = 1.0 if alpha is None else alpha

    z = np.empty((100, 1, 4), dtype=float)
    rgb = mcolors.colorConverter.to_rgb(fill_color)
    z[:,:,:3] = rgb
    z[:,:,-1] = (np.linspace(0., alpha, 100)[:,None])**(1.3)

    xmin, xmax, ymin, ymax = x.min(), x.max(), y.min(), y.max()
    im = ax.imshow(z, aspect='auto', extent=[xmin, xmax, ymin, ymax],
                   origin='lower', zorder=zorder)

    xy = np.column_stack([x, y])
    xy = np.vstack([[xmin, ymin], xy, [xmax, ymin], [xmin, ymin]])
    clip_path = Polygon(xy, facecolor='none', edgecolor='none', closed=True)
    ax.add_patch(clip_path)
    im.set_clip_path(clip_path)

    ax.autoscale(True)
    return line, im

def main():
    print("Hello from the plot program!")

    filenames = [
        "./../data/results/tcAsOfOmOneGapBos0.hdf5",
        "./../data/results/tcAsOfOmOneGapBos100.hdf5",
        "./../data/results/tcAsOfOmOneGapShortBos100.hdf5",
    ]
    plotTcAsOfOmega(filenames)
    exit()

    filename = "./../data/results/tcAsOfNphOneGapBos297.hdf5"
    plotPhaseDiagram(filename)


main()