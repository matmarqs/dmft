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
    filename = sys.argv[1]
    # loading data
    freq, ReG1, ImG1, ReG2, ImG2, ReG3, ImG3 = np.loadtxt(filename, unpack=True)
    spectral_func1 = spectral(ImG1)
    spectral_func2 = spectral(ImG2)
    spectral_func3 = spectral(ImG3)
    # plotting
    plt.plot(freq,  spectral_func1, label=r'$A_1(\omega)$')
    plt.plot(freq,  spectral_func2, label=r'$A_2(\omega)$')
    plt.plot(freq,  spectral_func3, label=r'$A_3(\omega)$')
    plt.xlim([-10.0, 10.0])
    plt.xlabel(r'$\omega$', fontsize=20)
    plt.ylabel(r'$A(\omega)$', fontsize=20)
    plt.legend(fontsize=16)
    plt.title(r'$A(\omega)$'+" of "+filename+" in the example")
    plt.savefig("specfun-examp.png", dpi=300, format='png', bbox_inches="tight")
    plt.clf()
    return 0


if __name__ == '__main__':
    main()
