#!/usr/bin/env python3

import numpy as np

from matplotlib import pyplot as plt
## comentar as 4 linhas abaixo caso nao tenha o LaTeX no matplotlib ###
from matplotlib import rc
plt.style.use('bmh')
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
#######################################################################

import sys

def spectral(ImG):
    return -ImG / np.pi

def main():
    file = sys.argv[1]
    # loading data
    freq, ReG, ImG = np.loadtxt(file, unpack=True)
    spectral_func = spectral(ImG)
    # plotting
    plt.plot(freq,  spectral_func, label=r'$A(\omega)$')
    #plt.xlim([-5.0, 5.0])
    plt.xlim([-5.0, 5.0])
    plt.xlabel(r'$\omega$', fontsize=20)
    plt.ylabel(r'$A(\omega)$', fontsize=20)
    plt.legend(fontsize=16)
    plt.title(r'$A(\omega)$'+" of "+file)
    plt.savefig("specfun.png", dpi=300, format='png', bbox_inches="tight")
    plt.clf()
    return 0


if __name__ == '__main__':
    main()
