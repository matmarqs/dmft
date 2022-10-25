#!/usr/bin/env python3

import os, subprocess
import numpy as np
from scipy.integrate import simpson as intg

dos_type = 'hypcub'  # 'bethe' or 'hypcub'
#t = 0.5 # D = 1 as energy unit, incomplete paper
t = 1.0     # t = 1 as energy unit, Bruno thesis
D = 2 * t
U = 1.7
#U = 1.75
#T = 0.0025
T = 1/7.2   # hypercubic, Bruno thesis
mu = U/2.0
max_err = 1e-4
alpha = 0.8
mesh1 = 'ppmesh1.dat'   # file for bath function mesh
mesh2 = 'ppmesh2.dat'   # file for pp self-energy mesh
nca = './nca-latest'
cix = 'cix=latest/cix.dat'
w0 = 0.01


def dos_bethe(E):   # density of states, bethe lattice
    return 2/(np.pi*D) * np.sqrt(1 - (np.clip(E, -D, D)/D)**2)
def dos_hypcub(E):     # density of states, hypercubic lattice
    return np.exp(-E*E/(t*t)) / (t * np.sqrt(np.pi))


def mix(new, old, alpha=0.5):
    return alpha * new + (1-alpha) * old


def broyden_init(v1, f1):
    v2 = v1 + alpha * f1
    df = []
    u = []
    return v2, df, u, 2


def broyden(v2, v1, f2, f1, df, u, m):
    if m <= 30:
        df_e = (f2 - f1) / np.linalg.norm(f2 - f1)
        df.append(df_e)
        df_arr = np.array(df, ndmin=2)
        dv = (v2 - v1) / np.linalg.norm(v2 - v1)
        u.append(alpha * df_e + dv)
        u_arr = np.array(u, ndmin=2)
        A = np.transpose(np.matmul(np.conj(df_arr), np.transpose(df_arr)))
        beta = np.linalg.inv(w0*w0 + A)
        c = np.conj(df_arr) @ f2
        v3 = v2 + alpha*f2 - c @ (beta @ u_arr)
    else:
        v3 = alpha * v2 + (1-alpha) * v1
    return v3, m+1

def main():
    if dos_type == 'bethe':
        dos = dos_bethe
    else:
        dos = dos_hypcub

    if not os.path.exists('out'):   # make sure the directory 'out' exists
        os.makedirs('out')

    # bath input file
    with open('delta.in', 'w') as delta_input:
        subprocess.run(['./gaumesh.py', mesh1], stdout=delta_input)
    freq = np.loadtxt('delta.in', unpack=True, ndmin=1)   # ndmin=2 to be a matrix
    ind = range(len(freq))
    # initial guess is Sigma = - U/20 * 1j
    Gf = np.array([intg(dos(freq) / (freq[i] + mu - U/20.0 * (-1j) - freq), freq) for i in ind])
    g0 = Gf     # g0 = [ 1/Gf + Sigma ]^{-1}, but Sigma=0 at initial guess
    #g0_1 = g0
    Delta = - np.imag(freq + mu - 1/g0) / np.pi
    np.savetxt('delta.in', np.transpose((freq, Delta)), fmt='%15.8e')
    # PP self-energies input file
    with open('sig.in', 'w') as sig_input:
        subprocess.run(['./gaumesh.py', mesh2], stdout=sig_input)
    #freq2 = np.loadtxt('sig.in', unpack=True, ndmin=1)   # ndmin=2 to be a matrix
    #freq2 = np.loadtxt('sig.in', unpack=True, ndmin=1)   # ndmin=2 to be a matrix
    #pp_sig = 0*freq2 + U/20.0
    #np.savetxt('sig.in', np.transpose((freq2, pp_sig, pp_sig, pp_sig)), fmt='%15.8e')

    # solve SIAM for the first time
    subprocess.call([nca, 'out=out', 'Sig=sig.in', 'Ac=delta.in', cix,
                    f'U={U}', f'T={T}', f'Ed={-mu}', 'max_diff=1e-6', 'max_steps=300'])
    Gimp = np.loadtxt('out/gloc.out', unpack=True, ndmin=2)   # ndmin=2 is to be a matrix
    g_imp = Gimp[1] + 1j * Gimp[2]
    Sig1 = 1/g0 - 1/g_imp

    Gf = np.array([intg(dos(freq) / (freq[i] + mu - Sig1[i] - freq), freq) for i in ind])
    g0 = 1/(1/Gf + Sig1)
    #fg1 = g0 - g0_1
    #g0, dg, ug, mg = broyden_init(g0_1, fg1)
    Delta = - np.imag(freq + mu - 1/g0) / np.pi
    np.savetxt('out/delta.loop', np.transpose((freq, Delta)), fmt='%15.8e')
    # solve SIAM for a second time
    subprocess.call([nca, 'out=out', 'Sig=out/Sigma.000', 'Ac=delta.in', cix,
                    f'U={U}', f'T={T}', f'Ed={-mu}', 'max_diff=1e-6', 'max_steps=300'])
    Gimp = np.loadtxt('out/gloc.out', unpack=True, ndmin=2)   # ndmin=2 is to be a matrix
    g_imp = Gimp[1] + 1j * Gimp[2]
    Sig2 = 1/g0 - 1/g_imp

    # v1 = Sig1
    f1 = Sig2 - Sig1
    Sig2, df, u, m = broyden_init(Sig1, f1)
    # v2 = Sig2

    # loop
    diff = np.max(np.abs(Sig2 - Sig1))
    while diff > max_err:
        Gf = np.array([intg(dos(freq) / (freq[i] + mu - Sig2[i] - freq), freq) for i in ind])
        #fg2 = 1/(1/Gf + Sig2) - g0
        #g0_aux = g0
        #g0, mg = broyden(g0, g0_1, fg2, fg1, dg, ug, mg)
        #g0_1 = g0_aux
        #fg1 = fg2
        Delta = - np.imag(freq + mu - 1/g0) / np.pi
        np.savetxt('out/delta.loop', np.transpose((freq, Delta)), fmt='%15.8e')
        # solve SIAM. 'Sig=out/Sigma.000' or 'Sig=sig.in'
        subprocess.call([nca, 'out=out', 'Sig=out/Sigma.000', 'Ac=out/delta.loop', cix,
                        f'U={U}', f'T={T}', f'Ed={-mu}', 'max_diff=1e-6', 'max_steps=300'])
        Gimp = np.loadtxt('out/gloc.out', unpack=True, ndmin=2)   # ndmin=2 is to be a matrix
        g_imp = Gimp[1] + 1j * Gimp[2]
        ######################################
        f2 = (1/g0 - 1/g_imp) - Sig2
        # mixing of self-energy
        Sig_aux = Sig2
        Sig2, m = broyden(Sig2, Sig1, f2, f1, df, u, m)
        Sig1 = Sig_aux
        f1 = f2
        diff = np.max(np.abs(Sig2 - Sig1))
        print(f'iteration {m}: diff = {diff}')


    # save the Green's function in the end
    np.savetxt('gf_dmft.out', np.transpose((freq, np.real(g_imp), np.imag(g_imp))),
                fmt='%15.8e')


if __name__ == '__main__':
    main()
