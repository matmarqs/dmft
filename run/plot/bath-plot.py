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

def main():
    file = sys.argv[1]
    # loading data
    freq, bath = np.loadtxt(file, unpack=True)
    # plotting
    plt.plot(freq,  bath, label=r'$\Delta(\omega)$')
    plt.xlim([-10.0, 10.0])
    plt.xlabel(r'$\omega$', fontsize=20)
    plt.ylabel(r'$\Delta(\omega)$', fontsize=20)
    plt.legend(fontsize=16)
    plt.title(r'$\Delta(\omega)$'+" of "+file)
    plt.savefig("bath.png", dpi=300, format='png', bbox_inches="tight")
    plt.clf()
    return 0


if __name__ == '__main__':
    main()
