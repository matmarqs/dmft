#!/usr/bin/env python3

import os, subprocess
import numpy as np
from gaumesh import genmesh

from matplotlib import pyplot as plt
## comentar as 4 linhas abaixo caso nao tenha o LaTeX no matplotlib ###
from matplotlib import rc
plt.style.use('bmh')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
#######################################################################

version = 'latest'  # version = '2004' or 'latest'
D = 1.0
V = 0.05
U = 0.27
Delta = np.pi * V**2 / (2 * D)
T = 1e-5
e0 = -U/2.0
#D = np.pi * V**2 / (2*Delta)    # half-bandwidth for calculations
b = 3.3 * D     # half-bandwith only for plots
nca = f'./nca-{version}'
cix = f'cix-{version}.dat'
figs_dir = f'figs_nca-{version}'
dx0  = [2e-7, 2e-5, 2e-4, 2e-2, 0.2]
fwhm = [2e-5, 2e-3, 2e-2, 2.0, 20.0]
x0   = [0.0,  0.0,  0.0,  0.0,  0.0]

def spectral(ImG):
    return -ImG / np.pi

def plot_specfun(file, figs_dir):
    # loading data
    data = np.loadtxt(file, unpack=True)
    specfun = spectral(data[2])
    # plotting
    plt.plot(data[0], specfun, label=r'$A(\omega)$')
    plt.xlim([-b, b])   # we plot for x in [-b, b]
    plt.xlabel(r'$\omega$', fontsize=20)
    plt.ylabel(r'$A(\omega)$', fontsize=20)
    plt.legend(fontsize=16)
    plt.title(r'$A(\omega)$'+" of "+file)
    plt.savefig(figs_dir+"/specfun.png", dpi=300, format='png', bbox_inches="tight")
    plt.clf()


def main():
    if not os.path.exists('out'):   # make sure the directory 'out' exists
        os.makedirs('out')

    mesh1 = np.array(genmesh(dx0, fwhm, x0, xmin=-2*b, xmax=2*b))
    mesh2 = np.array(genmesh(dx0, fwhm, x0, xmin=-4*b, xmax=4*b))
    delta = (np.abs(mesh1) <= D) * Delta
    np.savetxt('delta.in', np.transpose((mesh1, delta)), fmt='%15.8e')
    np.savetxt('sig.in', np.transpose((mesh2)), fmt='%15.8e')
    pp = mesh2 * 0 + Delta/2.0
    np.savetxt('sig.in', np.transpose((mesh2,pp,pp,pp)), fmt='%15.8e')

    subprocess.call([nca, 'out=out', 'Sig=sig.in', 'Ac=delta.in',
        f'cix={cix}', f'U={U}', f'T={T}', f'Ed={e0}',
        'max_diff=1e-8', 'max_steps=300'])

    if os.path.isfile('./cores.dat'):   # remove file ./cores.dat if it exists
        os.remove('./cores.dat')

    if not os.path.exists(figs_dir):   # make sure the directory figs_dir exists
        os.makedirs(figs_dir)
    plot_specfun('out/gloc.out', figs_dir)


if __name__ == '__main__':
    main()
