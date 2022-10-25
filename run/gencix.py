#!/usr/bin/env python3

# This script was originally in python2
# It was ported to python3 using the command '2to3 -w orig-gencix.py2'
# and then manual adjusts were made (they are marked "by Mateus")
# Also, the lines that generate 'ppspin.dat' were commented out
# It seems that this file was not needed.

########################################################
#                                                      #
# gencix.py                                            #
#                                                      #
# (C) by Kristjan Haule and David Jacob                #
#                                                      #
# previous atom_d.py scriptoriginally developed        #
# by Kristjan Haule for DMFT code                      #
# modified by David Jacob                              #
#                                                      #
# generates cix inputfile for OCA impurity solver      #
# from exact diagonalization of a multilevel impurity  #
# with general Coulomb interaction                     #
#                                                      #
########################################################

from scipy import *
from numpy import *     # by Mateus

import sys, re, os
import re,os
import copy
import getopt
import pickle

# wrapper to port Python2 to Python3, by Mateus
def cmp(a, b):
    return (a > b) - (a < b)

# wrapper to port Python2 to Python3, by Mateus
# cmp is an argument of list.sort() method
def cmp_to_key(mycmp):
    'Convert a cmp= function into a key= function'
    class K:
        def __init__(self, obj, *args):
            self.obj = obj
        def __lt__(self, other):
            return mycmp(self.obj, other.obj) < 0
        def __gt__(self, other):
            return mycmp(self.obj, other.obj) > 0
        def __eq__(self, other):
            return mycmp(self.obj, other.obj) == 0
        def __le__(self, other):
            return mycmp(self.obj, other.obj) <= 0
        def __ge__(self, other):
            return mycmp(self.obj, other.obj) >= 0
        def __ne__(self, other):
            return mycmp(self.obj, other.obj) != 0
    return K

def union(data1, data2):
    " Takes a union of two lists"
    res = data1
    for d in data2:
        if d not in res: res.append(d)
    return res

def overlap(data1, data2):
    " Checks if there is any overlap between data1 and data2"
    for d in data2:
        if d in data1:
            return True
    return False

def compress(groups):
    loopsDone = True
    while (loopsDone):
        loopsDone = False
        for i in range(len(groups)):
            if loopsDone: break
            for j in range(i+1,len(groups)):
                if loopsDone: break
                if overlap(groups[i],groups[j]):
                    groups[i] = union(groups[i],groups[j])
                    del groups[j]
                    loopsDone = True

    for g in groups: g.sort()

    groups.sort(key=cmp_to_key(lambda x,y: x[0]-y[0]))  # by Mateus

    return groups

def cprint(U):
    for i in range(shape(U)[0]):
        for j in range(shape(U)[1]):
            f = U[i,j]
            if abs(f)<1e-10: f = 0j
            print("%7.4f %7.4f*i " % (f.real, f.imag), end=' ')
        print()
    print()

def mprint(U):
    for i in range(shape(U)[0]):
        for j in range(shape(U)[1]):
            f = U[i,j].real
            if abs(f)<1e-10: f = 0
            print("%7.4f " % f, end=' ')
        print()
    print()

# Reads a generic Coulomb matrix U_ijkl from file
# Format:
# U_1111 U_1112 U_1113 ... U_1121 U_1122 ...
# U_1211 U_1212 U_1213 ...
# ...
# U_2111 U_2112 U_2113 ...
# ...
def ReadCoulombMatrix(nlevels,fname,UC):
    f = open(fname, 'r')
    for i in range(nlevels):
        for j in range(nlevels):
            line = f.readline()
            if(len(line.split())<nlevels*nlevels): break
            KL=0
            for k in range(nlevels):
                for l in range (nlevels):
                    UC[i,j,k,l] = float(line.split()[KL])
                    KL+=1
    f.close()

# Read generic Coulomb matrix U_ijkl in sparse matrix Form
# Each line contains one matrix element
# Format: i j k l U_ijkl
def ReadCoulombMatrixSparse(nlevels,fname,UC):
    f = open(fname, 'r')
    line = f.readline()
    while f:
        if(len(line.split())<5): break
        i = int(line.split()[0])-1
        j = int(line.split()[1])-1
        k = int(line.split()[2])-1
        l = int(line.split()[3])-1
        UC[i,j,k,l] = float(line.split()[4])
        line = f.readline()
    f.close()

# General many-body operator
# implementing one and two-body interactions
# acting on Slater determinants
class operateLS(object):
    def __init__ (self, Nband):
        self.Nband = Nband
        self.baths = 2*Nband

        self.N = self.baths

        self.mask=[]
        for i in range(self.N): self.mask.append(1<<i);

        self.bi=[] # band index
        self.sz=[] # sz
        for i in range(self.Nband):
            self.sz.append(1); self.sz.append(-1)
            self.bi.append(i); self.bi.append(i)

        self.mask_d = []
        self.mask_u = []
        for i in range(self.Nband):
            self.mask_u.append(self.mask[2*i])
            self.mask_d.append(self.mask[2*i+1])

    def printn(self, state):
        sstate=''
        for i in range(self.Nband):
            if (state & self.mask_u[i]) and (state & self.mask_d[i]) : sstate += '2'
            elif (state & self.mask_u[i]) : sstate += 'u'
            elif (state & self.mask_d[i]) : sstate += 'd'
            else : sstate += '0'
            #sstate += ' '
        return sstate

    def Nel(self, state):
        n=0
        for k in self.mask:
            if (k&state): n+=1
        return n

    def Sz(self, state):
        nu = 0
        nd = 0
        for i in range(self.Nband):
            if state&self.mask_u[i] : nu += 1
            if state&self.mask_d[i] : nd += 1
        return nu-nd

    def occup(self, state):
        """ gives a list of occupancies per band [n_{band1},n_{band2},...]
        """
        oc=[]
        for i in range(self.Nband):
            ne = 0
            if state & self.mask_u[i] : ne += 1
            if state & self.mask_d[i] : ne += 1
            oc.append(ne)
        return oc

    def sp_occup(self, state):
        """ Spin-resolved occupation numbers
        gives a list of occupancies per bath (i.e. spin band)  [n_{band1},n_{band2},...]
        """
        oc=[]
        for i in range(self.baths):
            ne = 0
            if state & self.mask[i] : ne = 1
            oc.append(ne)
        return oc

    def sign(self, state, mask_min, mask_max):
        """ Sign when electron hops from mask_min to mask_max
        Counts number of electrons between the two spaces
        """
        # mask will run between mask_min to mask_max
        mask = mask_min<<1
        n=0                       # number of electrons between mask_min and mask_max
        while (mask<mask_max):    # loop to mask_max
            if (mask&state): n+=1 # found electron between the two places
            mask = mask<<1        # increment the mask
        return 1-2*(n%2)          # (-1)^n

    def sign_(self, state, mask_max):
        """ Sign when electron is added to the state (from the left)
        """
        # mask will run between mask_min to mask_max
        mask = 1
        n=0           # number of electrons between mask_min and mask_max
        while (mask<mask_max):    # loop to mask_max
            if (mask&state): n+=1 # found electron between the two places
            mask = mask<<1        # increment the mask
        return 1-2*(n%2)          # (-1)^n

    def N_el_before(self, state, i):
        n=0
        for q in range(i):
            if self.mask[q]&state: n+=1
        return n

    def S2(self, state):
        l2p1 = self.Nband
        sts=[]
        # diagonal part
        dd=0;
        for ilz in range(l2p1):
            up=0; dn=0
            if self.mask_u[ilz] & state: up = 1
            if self.mask_d[ilz] & state: dn = 1
            # if only up or only down in certain lz
            if up+dn==1: dd += 0.5
        # Sz^2
        fct = (0.5*self.Sz(state))**2 + dd
        # store diagonal
        sts.append([state,fct])
        # off diagonal
        for ilz in range(l2p1):
            im1 = self.mask_u[ilz]
            im2 = self.mask_d[ilz]
            ib1 = bool(state & im1)
            ib2 = bool(state & im2)
            if ib1 and not ib2: # S^-_i gives nonzero
                isig = self.sign(state, min(im1,im2), max(im1,im2))
                istate = state^im1^im2
                for jlz in range(l2p1):
                    if (ilz==jlz): continue
                    jm1 = self.mask_d[jlz]
                    jm2 = self.mask_u[jlz]
                    jb1 = bool(state & jm1)
                    jb2 = bool(state & jm2)
                    if jb1 and not jb2: # S^+_j gives nonzero
                        jsig = self.sign(istate, min(jm1,jm2), max(jm1,jm2))
                        jstate = istate^jm1^jm2
                        sts.append([jstate, isig*jsig])
        return sts

    def Fp(self, state, b, s):
        """ This implements psi^dagger_{b,s} operator acting on state
        indexes are:
          b-band
          s-spin
        """
        if (s==0): mask = self.mask_d[b]
        elif (s==1): mask = self.mask_u[b]
        else:
            print('spin should be 0 or 1!')
            sys.exit(1)

        if state&mask: return (0,1)  # This state is already occupied
        newstate = state^mask
        sig = self.sign_(state, mask)
        return (newstate, sig)

    # Coulomb interaction given by Coulomb matrix UC[i,j,k,l]
    # acting on many-body state (Slater determinant) 'state'
    # gives back an array of slater determinants that have non-zero matrix elements with the original 'state'
    # and the corresponding matrix element
    def CoulombU(self, state, UC ):
        sts=[]
        ni=-1
        for i in range(self.baths):
            if not(self.mask[i]&state) : continue
            ni+=1
            #ni = self.N_el_before(state,i)
            state1 = state^self.mask[i]
            m1 = self.bi[i]
            s1 = self.sz[i]
            nj=-1
            for j in range(self.baths):
                if not(self.mask[j]&state1) : continue
                nj+=1
                #nj = self.N_el_before(state1,j)
                # here we have: mask[i]&state && mask[j]&state
                state2 = state1^self.mask[j]
                m2 = self.bi[j]
                s2 = self.sz[j]
                for a in range(self.baths):
                    if self.mask[a]&state2 or self.sz[a]!=s2 : continue # sz_a == sz_j
                    na = self.N_el_before(state2,a)
                    state3 = state2^self.mask[a]
                    m3 = self.bi[a]
                    s3 = self.sz[a]
                    for b in range(self.baths):
                        if self.mask[b]&state3 or self.sz[b]!=s1: continue # sz_b == sz_i
                        nb = self.N_el_before(state3,b)
                        state4 = state3^self.mask[b]
                        m4 = self.bi[b]
                        s4 = self.sz[b]
                        #sign = 1-2*((ni+nj+na+nb)%2)
                        sign = (-1)**(ni+nj+na+nb)
                        U0 = sign*UC[m4,m3,m1,m2]
                        if (abs(U0)>1e-6): sts.append([state4, U0, m4, m3, m1, m2])
        return sts


def baseN(Nband, prop):
    Ntot = len(prop)
    wstates=[]

    for n1 in range(Nband*2+1):
        for sz1 in range(-n1,n1+1,2):
            states=[]
            for i in range(Ntot):
                if prop[i][0]==n1 and prop[i][1]==sz1:
                    states.append(i)
            if (len(states)>0): wstates.append([n1, sz1, states])
    return wstates


def list_to_string(x):
    return str(array(x).flatten().tolist())


def analizeGroups(A, small = 1e-6):
    groups=[]
    for i in range(shape(A)[0]):
        nonz=[]
        for j in range(shape(A)[1]):
            if abs(A[i,j])>small : nonz.append(j)
        if (len(nonz)>0): groups.append(nonz)

    groups0 = compress(groups)

    groups=[]
    for i in range(shape(A)[1]):
        nonz=[]
        for j in range(shape(A)[0]):
            if abs(A[j,i])>small : nonz.append(j)
        if (len(nonz)>0): groups.append(nonz)

    groups1 = compress(groups)
    return (groups1,groups0)

def coupled(A, groups0, groups1, small = 1e-6):
    #ng0 = len(array(groups0).flatten().tolist())
    #ng1 = len(array(groups1).flatten().tolist())

    fpair = [-1]*len(groups0)
    #pairs=[]
    for ii,ig0 in enumerate(groups0):
        nonz=[]
        for ir0 in ig0:
            for q in range(shape(A)[1]):
                if abs(A[ir0,q])>small : nonz.append(q)
        for jj,jg1 in enumerate(groups1):
            if overlap(nonz,jg1):
                #pairs.append([ii,jj])
                fpair[ii] = jj
    return fpair

def comp(x, y):
    if x[2]!=y[2]: return int(x[2]-y[2])
    else:
        if (x[3]==y[3]): return 0
        elif (x[3]<y[3]): return -1
        else: return 1

def Diagonalize(Ham, small=1e-6):
    """ Diagonalization is done in blocks. This is not because of efficiency but because
    the resulting eigenvectors must not mix states of direct base if not absolutely necessary.
    If brute force diagonalization is used in large scale problems, eigenvectors can be seriously
    mix direct states with different symmetry.
    """
    # Check block structure of Hamiltonian
    # States which are mixed in Ham will have the same blck[i]
    ndim = len(Ham)
    blck=list(range(ndim))
    for i in range(ndim):
        for j in range(i+1,ndim):
            if (abs(Ham[i][j])>small): blck[j]=blck[i]
    # Having blck[i] a new array block[:][:] is creates, which contains indexes to all blocks
    # for example [[1,2,3],[4,5,6]] for Hamiltonian containing two blocks
    block=[]
    for i in range(ndim):
        bb=[]
        for j in range(ndim):
            if blck[j]==i: bb.append(j)
        if len(bb)>0:
            block.append(bb)

    #print 'blocks=', block

    # Here we go over all blocks and diagonalize each one.
    eigv=[] # contains all eigenvalues
    eigx=[] # contains full eigenvectors
    for bl in block:
        hs = zeros((len(bl),len(bl)), dtype=complex)
        for i,ib in enumerate(bl):
            for j,jb in enumerate(bl):
                hs[i,j] = Ham[ib,jb]

        eig1 = linalg.eigh(hs)

        # Checking if it is comple!
        for l in range(len(eig1[1])):
            ime = sum([abs(eig1[1][u,l].imag) for u in range(len(eig1[1]))])
            if ime>1e-7: print('TROUBLES!!! Complex eigenvector! You sould improve that!')

        # Creating a big eigenvector with all components
        for l in range(len(eig1[1])):
            large_eig=zeros(ndim, dtype=complex)
            small_eig = eig1[1][:,l]
            for m,mb in enumerate(bl):  large_eig[mb] = small_eig[m]
            eigx.append(large_eig)
        eigv += eig1[0].tolist()

    # Now we need to sort eigenvectors and eigenvalues
    # index is created for sorting
    indx=list(range(ndim))
    indx.sort(key=cmp_to_key(lambda a,b: cmp(eigv[a],eigv[b])))     # by Mateus
    # and actual sorting is performed
    seigv=[]
    seigx=[]
    for i in range(ndim):
        seigv.append(eigv[indx[i]])
        seigx.append(eigx[indx[i]])

    # Eigenvectors should be in the form Ham*v[:,i] = w[i]*v[:,i]
    # which means that we need to transpose the list of eigenvectors
    seigx = array(seigx).transpose()

    # We also do a brute force diagonalization just to check if something goes wrong with block diagonalization
    # Note that the two resulting eigensystems are not necessary the same due to freedom in choosing eigenvectors
    #eig = symeig.symeig(Ham)

    # If eigenvalues from block diagonalization and full diagonalization are different, something is wrong
    #if sum(map(abs,eig[0]-array(seigv)))>small:  print 'TEZAVE!'

    return [seigv,seigx]

def EquivalentBaths(Eimp):
    """ Finds which baths are equivalent from impurity levels"""
    wE = [(i,Eimp[i]) for i in range(len(Eimp))]
    #print 'wE=', wE
    kbths=[]
    while len(wE)>0:
        En = wE[0][1]
        j=0
        rr=[]
        while j < len(wE):
            if abs(En-wE[j][1])<1e-10:
                rr.append(wE[j][0])
                del wE[j]
            else: j+=1
            #print 'w', j, rr, En, wE[j][1]
        kbths.append(rr)

    bathis=list(range(len(Eimp)))
    for ik,k in enumerate(kbths):
        for ij in k: bathis[ij]=ik

    #print 'kbths=', kbths
    Ed=[]
    for ik,k in enumerate(kbths):
        Ed.append(Eimp[k[0]])

    return (bathis,kbths,Ed)

def thesame(mx,my,small=1e-3):
    if list(mx.keys()) != list(my.keys()): return False
    for k in list(mx.keys()):
        if abs(mx[k]-my[k])>small: return False
    return True

def VEquivalentStates(mps,ind):
    """ Finds which states have the same bubbles """
    wx = [(i,mps[i]) for i in range(len(mps))]
    iequiv=[]
    while len(wx)>0:
        mx = wx[0][1]
        j=0
        rr=[]
        while j < len(wx):
            if thesame(mx,wx[j][1]):
                rr.append(wx[j][0])
                del wx[j]
            else: j+=1
        iequiv.append(rr)

    for ik in range(len(iequiv)):
        for ij in range(len(iequiv[ik])):
            iequiv[ik][ij] = ind[iequiv[ik][ij]]

    return iequiv

def AverageBubbles(tmps):
    """ Compute average over 'almost' equivalent states """
    trmp=[]
    for mps in tmps:
        all_keys=[]
        for mp in mps:
            all_keys = union(all_keys,list(mp.keys()))

        rmp={}
        for k in all_keys:
            sm=0.0
            for mp in mps:
                if k in mp:
                    sm += mp[k]
            #sm/=len(mps)
            rmp[k]=sm
        trmp.append(rmp)
    return trmp


def EquivalentStates(ipE, ipN):
    iequiv=[]
    equiv = list(range(len(ipE)))
    leq=0
    ju=0
    Nmax = ipN[-1]
    for Ni in range(Nmax+1):
        # all states of the same N are in the interval [ju,je]
        je=ju
        while je<len(ipN) and ipN[je]==Ni: je+=1

        ind = list(range(ju,je))
        ind.sort(key=cmp_to_key(lambda x,y: cmp(ipE[x],ipE[y])))    # by Mateus

        #print Ni
        i0=0
        while (i0<len(ind)):
            Ec = ipE[ind[i0]]
            ieq=[]
            #print 'Ec=', Ec
            while i0<len(ind) and abs(ipE[ind[i0]]-Ec)<1e-10:
                #print ind[i0], ipE[ind[i0]], leq
                equiv[ind[i0]] = leq
                ieq.append(ind[i0])
                i0+=1
            leq += 1
            iequiv.append(ieq)
            #print
        #print
        ju=je
    return (equiv, iequiv)

def RenumberStates(pseudostates, Enes, wstates, S2ws):
    # renumbers states such that each of 1024 states has unique index
    # also remembers energy and N for each state
    ij=0
    puniq={}
    ipuniq=[]
    ipE=[]
    ipN=[]
    ipS=[]
    for ii,iwp in enumerate(pseudostates):
        wdim = len(Enes[ii])
        for j in range(wdim):
            puniq[(ii,j)]=ij
            ipuniq.append((ii,j))
            ipE.append(Enes[ii][j])
            ipS.append(S2ws[ii][j])
            wN = sum(wstates[iwp[0]][0])
            ipN.append(wN)
            ij+=1
    return (puniq, ipE, ipN, ipS)


def CreateEmpty3D_Dict(n0,n1,n2):
    return [[[{} for i2 in range(n2)] for i1 in range(n1)] for i0 in range(n0)]
def CreateEmpty2D_Dict(n0,n1):
    return [[{} for i1 in range(n1)] for i0 in range(n0)]

#
# Main program starts here...
#

if __name__ == '__main__':
    """ Help here"""

    nlevels=1    # number of impurity levels
    n=[0,1,2]    # occupanices used for OCA
    Ncentral=[1] # OCA diagrams are selected such that central occupancy is in Ncentral
    qOCA=1       # OCA diagrams are computed
    Eoca=10.     # Energy window for OCA diagrams
    mOCA=1e-3    # matrix element for OCA should be greater than that
    Ewindow = [-1000,1000]
    max_M_size=6
    Ucfile = 'none' # external file for reading in generic Coulomb matrix U_ijkl

    # Parameters for model interaction used when
    # generic interaction matrix Uc is not specified (default)
    U0 = 2.0   # direct repulsion between electrons in same level
    U1 = 1.0   # direct repulsion between electrons in different levels
    JH = 0.5   # Hund's rule coupling

    # Array of impurity levels
    Eimp = [-1.0]

    Nitt= 1 # To figure out the symmetry, we iterate Nitt times

    args = sys.argv[1:]
    if ('-h' in args) or ('--help' in args):
        print("""Code for generating impurity cix file for a multilevel impurity with generic interaction
                 The output cix-files can be used for OCA or CTQMC solvers
                 Output files:
                    out.cix       -- the oca cix file
                    impurity.cix  -- the ctqmc cix file
                 The input is:
                    nlevels       -- number of impurity levels
                    Eimp          -- list of impurity levels (i.e., [0,0,0...] )
                    Ucfile        -- file name of file with generic Coulomb interaction matrix U_ijkl
                    U0            -- direct Coulomb repulsion for electrons in same level (intra-level repulsion)
                    U1            -- direct Coulomb repulsion for electrons in different levels (inter-level repulsion)
                    JH            -- Hund's rule coupling
                    n             -- list of occupancies to be taken into account
                    NCentral      -- Central occupancy for OCA diagrams
                    qOCA          -- OCA diagrams included if qOCA=1
                    Eoca          -- OCA diagrams cutoff: if any of the atomic states has for Eoca higher energy than the ground state for particular occupancy, the diagram is dropped
                    mOCA          -- If matrix element is smaller than mOCA, the diagram is dropped
                    Ewindow       -- Energy window for the states kept (used in ctqmc only)
                    max_M_size    -- maximum matrix size kept for ctqmc
                 """)
        sys.exit(0)
    for arg in args:
        if os.path.isfile(arg):
            exec(compile(open(arg, "rb").read(), arg, 'exec'))
            print('Executed file', arg)
        else:
            exec(arg)

    # Read generic Coulomb interaction matrix from file
    # and build full Coulomb interaction matrix Uc
    if( Ucfile != 'none' ):
        Uc = zeros((nlevels,nlevels,nlevels,nlevels), dtype=float)
        ReadCoulombMatrixSparse(nlevels,Ucfile,Uc)

    print('nlevels = ', nlevels)
    print('Eimp = ', '[','%f, '*len(Eimp) % tuple(Eimp),']')
    print('n = ', n)
    if( Ucfile == 'none' ):
        print('Model interaction parameters:')
        print('U0 = ', U0)
        print('U1 = ', U1)
        print('JH  = ', JH)
    else:
        print('Coulomb interaction matrix read from file', Ucfile)
        print('Uc = ')
        print(Uc)
    print()
    print('qOCA = ', qOCA)
    print('Eoca = ', Eoca)
    print('mOCA = ', mOCA)
    print('Ewindow = ', Ewindow)
    print('max_M_size=', max_M_size)

    # generate one electron basis of impurity (bath = impurity spin level)
    baths=[]
    for b in range(nlevels):
        for s in range(2):
            baths.append([b,s])

    (bathi, kbth, Ed) = EquivalentBaths( [Eimp[b[0]] for b in baths] )

    Ntot = 2**(len(baths)) # size of the direct base

    op = operateLS(nlevels) # initialize operator class

    # some properties of integers which will serve as a direct base - partial occupancy and Sz
    prop=[]
    for i in range(Ntot):
        occ = op.occup(i)
        prop.append([sum(occ), op.Sz(i),occ])
    # creates direct base from integers having correct properties
    # wstates contains [N, Sz, [all direct states with this N and Sz]]
    wstates = baseN(nlevels,prop)

    indx={}
    for ni,ns in enumerate(wstates):
        indx[(ns[0],ns[1])] = ni  # index according to N and Sz of the state
        print(ns[0], ns[1], ni)

    kindx = list(indx.keys())

    print('Stage0: Exact diagonalization of the atom', file=sys.stderr)

    Ene=[] # Energy
    Te=[]  # eigenvectors
    S2w=[] # Spin
    for ni,ns in enumerate(wstates):
        print('n=', ns[0], 'sz=', ns[1], 'dim=', len(ns[2]))
        states = ns[2]

        # Computes matrix of S^2
        S2 = zeros((len(states),len(states)),dtype=complex)
        for js,st in enumerate(states):
            stn = op.S2(st)
            for ps in stn:
                ii = ps[0]
                iu = states.index(ii)
                S2[js,iu] += ps[1]

        # Now we compute the Hamiltonian in direct basis
        Ham = zeros((len(states),len(states)),dtype=complex)

        for js,st in enumerate(states):
            #
            # on-site energies contain crystal-field splittings
            #
            occ = op.occup(st)
            for ic,oc in enumerate(occ):
                if (js>=len(states)): print('Tezave!')
                Ham[js,js] += oc*Eimp[ic]
            #
            # Coulomb inbteraction in direct basis
            #
            # 1. generic Coulomb interaction U_ijkl
            if Ucfile != 'none':
                # all direct states cst generated
                # by acting with Coulomb on direct state st
                cst = op.CoulombU(st, Uc )
                for cs in cst:
                    ii = cs[0]
                    U0 = cs[1]
                    iu = states.index(ii)
                    Ham[js,iu] += 0.5*U0
            # 2. Model Coulomb interaction
            # containing only direct repulsion
            # and Hunds rule coupling
            else:
                occ = op.occup(st)
                ndd = 0 # number of doubly occupied levels
                # Compute direct repulsion (diagonal in direct basis!)
                for i1,oc1 in enumerate(occ):
                    if oc1==2: ndd += 1
                    for i2,oc2 in enumerate(occ):
                        # intra-level repulsion
                        if i1==i2: Ham[js,js] += 0.5 * U0 * ( oc1 * oc2 - oc1 )
                        # inter-level repulsion
                        if i1!=i2: Ham[js,js] += 0.5 * U1 * oc1 * oc2
                    Ham[js,js]
                # Compute Hunds rule coupling
                # -JH*S^2 + JH*nd - 1/4*J*nd^2 - JH*ndd
                # nd  = total d-shell occupation,
                # ndd = number of doubly occupied d-levels
                cst = op.S2(st)
                for cs in cst:
                    ii = cs[0]
                    iu = states.index(ii)
                    Ham[js,iu] += -JH*cs[1]
                nd = sum(occ)
                # Kristjan's original code:
                #Ham[js,js] += JH*nd*(1-0.25*nd)
                # not quite correct: one has to also substract number of doubly occupied levels
                # Hence corrected code:
                Ham[js,js] += JH*(nd*(1-0.25*nd)-ndd)

        # Diagonalization of Hamiltonian
        eig = Diagonalize(Ham)  # Block diagonalization
        T0 = eig[1]             # eigen vectors of Hamiltonian

        # Here we compute matrix of S^2 in eigen basis.
        # Should be diagonal without spin-orbit coupling
        S2e = matrix(conj(T0.transpose())) * S2 * matrix(T0)

        # Check if S^2 is diagonal
        printS=False
        for i0 in range(shape(S2e)[0]):
            for i1 in range(shape(S2e)[1]):
                if i0!=i1 and abs(S2e[i0,i1])>1e-6 :
                    print('Warning: S^2 not diagonal! ', i0, i1, S2e[i0,i1])
                    printS=True
        if printS:
            mprint(S2e)
            print()
            print(eig[0])
            print()
            mprint(Ham)

        Ene.append(eig[0])
        Te.append(eig[1])
        # Spin is computed using formula s(s+1)
        S2w.append([0.5*int(round(-1+sqrt(1+4*S2e[i,i].real))) for i in range(len(S2e))])

    # Here we create index for psi^dagger
    iFi = zeros((len(wstates),len(baths)),dtype=int)
    for ni,ns in enumerate(wstates):
        for ib,be in enumerate(baths):
            st = (ns[0]+1, ns[1] + 2*be[1]-1)  # (n+1,sz+s)
            if st in kindx:
                iFi[ni,ib] = indx[st]
    wgr=[]
    for iw in range(len(wstates)): wgr.append([])

    print('Stage1: Computing F^ in direct base', file=sys.stderr)
    # Below we compute matrix elements of F^ (FKP)
    kindx = list(indx.keys())
    FKP = []
    for ni,ns in enumerate(wstates):
        states = ns[2]
        bFp=[]

        for ib,wib in enumerate(baths):
            b = wib[0]
            s = wib[1]

            inew = iFi[ni,ib]

            if inew==0:
                bFp.append([])
                continue

            #print 'here I am...'

            newstates = wstates[inew][2]

            Fp = zeros((len(states), len(newstates)), dtype=float)

            for js,st in enumerate(states):

                (newst, sig) = op.Fp(st, b, s)

                if newst>0:
                    ii = newstates.index(newst)
                    Fp[js,ii] += sig

            Fpn = matrix(conj(Te[ni].transpose())) * matrix(Fp) * matrix(Te[inew])

            # Set to zero small values
            for i0 in range(shape(Fpn)[0]):
                for i1 in range(shape(Fpn)[1]):
                    if abs(Fpn[i0,i1])<1e-4: Fpn[i0,i1]=0.0

            gr = analizeGroups(Fpn)
            # |inew> = F^+ |ni>
            wgr[ni] += gr[0]   # which states are coupled by F^ in |ni>
            wgr[inew] += gr[1] # which states are coupled by F^ in |inew>

            bFp.append(Fpn)

        FKP.append(bFp)

    print('Stage2: Compressing F^+ according to its block diagonal form', file=sys.stderr)

    for i in range(len(wstates)): wgr[i] = compress(wgr[i])

    print('Stage3: Renumbers states -- creates superstates for ctqmc', file=sys.stderr)

    # Here we store ground state and N for each superstates to be used later for sorting
    tstates=[]
    for iw in range(len(wstates)):
        Nt = sum(wstates[iw][0])
        for ip in range(len(wgr[iw])):
            Eg = Ene[iw][wgr[iw][ip][0]]
            tstates.append([iw,ip,Nt,Eg])
    tstates.sort(key=cmp_to_key(comp))  # by Mateus
    # tstates now contains [index-to-wstates, index-to-state-inside-wstates, N, E]

    # superstates == pseudostates are defined
    pseudostates=[]
    indpseudo={}
    jj=0
    for st in tstates:
        iw = st[0]
        ip = st[1]
        pseudostates.append([iw,ip])
        indpseudo[(iw,ip)] = jj
        jj+=1


    iFi_inside=[]
    for iw in range(len(wstates)):
        biFi=[]
        for ib in range(len(baths)):
            ifi = iFi[iw,ib]
            if ifi>0:
                fpair = coupled(FKP[iw][ib], wgr[iw], wgr[ifi])
                biFi.append(fpair)
            else:
                biFi.append([])
        iFi_inside.append(biFi)

    #print 'pseu=', pseudostates

    # creates arrays containing Energy, occupancy and index table for all superstates
    iFinal = zeros((len(pseudostates),len(baths)),dtype=int)
    Enes = []
    S2ws = []
    Occ=[]
    for ii,iwp in enumerate(pseudostates):
        iw = iwp[0]
        ip = iwp[1]
        wstate = wstates[iw]
        group = wgr[iw][ip]

        for ib in range(len(baths)):
            ifi = iFi[iw,ib]
            ifinal = -1
            if (ifi>0):
                ifi_ins = iFi_inside[iw][ib][ip]
                if ifi_ins>=0:
                    ifinal = indpseudo[(ifi,ifi_ins)]
            iFinal[ii,ib] = ifinal

        Ens=[]
        occ=[]
        S2s=[]
        for iq,q in enumerate(group):
            Ens.append(Ene[iw][q])
            occ.append(wstate[0])
            S2s.append(S2w[iw][q])

        Enes.append(Ens)
        Occ.append(occ)
        S2ws.append(S2s)

    print('Stage4: F^dagger matrices between superstates evaluated', file=sys.stderr)
    # creates small F^dagger matrices between superstates
    maxs = 0
    rFKP = []
    for ii,iwp in enumerate(pseudostates):
        iw = iwp[0]
        ip = iwp[1]
        wstate = wstates[iw]
        group = wgr[iw][ip]

        bFKP=[]
        for ib in range(len(baths)):
            ifi = iFi[iw,ib]
            ifinal = iFinal[ii,ib]
            if (ifi>0): ifi_ins = iFi_inside[iw][ib][ip]

            if ifinal>=0:
                M = zeros((len(group),len(wgr[ifi][ifi_ins])),dtype=float)
                for ii0,i0 in enumerate(group):
                    for jj0,j0 in enumerate(wgr[ifi][ifi_ins]):
                        M[ii0,jj0] = FKP[iw][ib][i0,j0].real
                        if abs(FKP[iw][ib][i0,j0].imag)>1e-10:
                            print('TROUBLE!!! Complex eigenvector!!!', file=sys.stderr)
                if max(shape(M)) > maxs : maxs = max(shape(M))
            else:
                M = array([])
            bFKP.append(M)
        rFKP.append(bFKP)

    ############################################################################################
    ######## The part of the code between this symbols,  generates input for OCA solver  #######
    ############################################################################################
    print('Stage5: Renumber states for oca - each atomic states has unique number', file=sys.stderr)
    # renumbers states such that each of 1024 states has unique index
    # also remembers energy and N for each state
    (puniq, ipE, ipN, ipS) = RenumberStates(pseudostates, Enes, wstates, S2ws)

    # bubbles will be accessed by bubbles[ib][ii][ij]
    # where ib is integer running over all baths, ii is integer running over all possible states
    # and ij is integer which can be accessed by F^{+,ib}|ii>
    bubbles = CreateEmpty2D_Dict(len(kbth), len(ipE))
    FpF = zeros((len(ipE),len(kbth)))
    smallb=1e-4
    for ib,bs in enumerate(kbth): # over all different baths
        bubl=bubbles[ib]
        for ii,iwp in enumerate(pseudostates): # over all pseudostates
            for iib in bs: # over equivalent baths
                ifinal = iFinal[ii,iib]
                if ifinal>=0:
                    dims = shape(rFKP[ii][iib])
                    for i0 in range(dims[0]):
                        istr = puniq[(ii,i0)]
                        for j0 in range(dims[1]):
                            iend = puniq[(ifinal,j0)]
                            FpF[istr][ib] += rFKP[ii][iib][i0,j0]**2
                            if (abs(rFKP[ii][iib][i0,j0])>smallb):
                                if iend in bubl[istr]:
                                    bubl[istr][iend] += rFKP[ii][iib][i0,j0]**2
                                else:
                                    bubl[istr][iend] = rFKP[ii][iib][i0,j0]**2



    (equiv, iequiv) = EquivalentStates(ipE, ipN)
    # Now we have all the bubbles.
    # We need to find which states are equivalent
    for tt in range(Nitt):

        print(tt, '. iteration')

        # We want to merge bubbles that we believe are equivalent
        ebubbles=[]
        for ib,bs in enumerate(kbth):
            bubl = bubbles[ib]

            nbubl=[]
            for i in range(len(bubl)): nbubl.append({})

            for i0 in range(len(bubl)):
                for i1 in list(bubl[i0].keys()):
                    if equiv[i1] in nbubl[i0]:
                        nbubl[i0][equiv[i1]] += bubl[i0][i1]
                    else:
                        nbubl[i0][equiv[i1]] = bubl[i0][i1]
            ebubbles.append(nbubl)

        # Here we check if the states, which we idenfified above as equivalent,
        # really are equivalent, i.e., have the same type of bubble
        new_iequiv=[]
        back_bubs=[]
        for ii in iequiv:
            cbubs=[]
            for ij in ii:
                cbub={}
                for ib in range(len(kbth)):
                    cbub.update(ebubbles[ib][ij])
                cbubs.append(cbub)

            abubs = AverageBubbles([[ebubbles[ib][ij] for ij in ii] for ib in range(len(kbth))])

            ieqs = VEquivalentStates(cbubs,ii)

            back_bubs.append(abubs)

            new_iequiv += ieqs


        new_equiv=list(range(len(equiv)))
        for i,ii in enumerate(new_iequiv):
            for j in ii: new_equiv[j]=i

        if len(iequiv)==len(new_iequiv) or tt+1==Nitt: break
        equiv = new_equiv
        iequiv = new_iequiv


    if qOCA:  # Here we add second order OCA diagrams to NCA bubbles
        # Second order diagramms
        # OCAdiag will be accessed by OCAdiag[ib1][ib2][ii][ij]
        # where ib1, ib2 is integer running over all baths, ii is integer running over all possible states
        # and ij is integer which can be accessed by F^{+,ib}|ii>
        OCAdiag = CreateEmpty3D_Dict(len(kbth), len(kbth), len(iequiv))

        for ib1,bs1 in enumerate(kbth):     # over all baths ones
            for ib2,bs2 in enumerate(kbth): # over all baths twice
                OCAs=OCAdiag[ib1][ib2]
                for ii,iwp in enumerate(pseudostates): # over all pseudostates
                    for iib1 in bs1:        # over equivalent baths ones
                        for iib2 in bs2:    # over equivalent baths twice
                            ifinal_j = iFinal[ii,iib1]
                            ifinal_l = iFinal[ii,iib2]
                            if ifinal_j>=0 and ifinal_l>=0:
                                ifinal_k = iFinal[ifinal_j,iib2]
                             # print 'Estou aqui ifinal_l , iib1 = ', ifinal_l , iib1 , iFinal[1,iib1]
                                ifinal_k2 = iFinal[ifinal_l,iib1]
                                if ifinal_k>=0 and ifinal_k2>=0 and ifinal_k==ifinal_k2:
                                    Fij = rFKP[ii][iib1]
                                    Fil = rFKP[ii][iib2]
                                    Fjk = rFKP[ifinal_j][iib2]
                                    Flk = rFKP[ifinal_l][iib1]
                                    (dims_i, dims_j)  = shape(Fij)
                                    (dims_i2,dims_l)  = shape(Fil)
                                    (dims_j2, dims_k) = shape(Fjk)
                                    (dims_l2,dims_k2) = shape(Flk)
                                    if dims_i != dims_i2: print('Problem i')
                                    if dims_j != dims_j2: print('Problem j')
                                    if dims_l != dims_l2: print('Problem l')
                                    if dims_k != dims_k2: print('Problem k')

                                    for i0 in range(dims_i):
                                        iu = equiv[puniq[(ii,i0)]]
                                        for j0 in range(dims_j):
                                            ju = equiv[puniq[(ifinal_j,j0)]]
                                            if (abs(Fij[i0,j0])<smallb): continue
                                            for k0 in range(dims_k):
                                                ku = equiv[puniq[(ifinal_k,k0)]]
                                                if (abs(Fjk[j0,k0])<smallb): continue
                                                for l0 in range(dims_l):
                                                    lu = equiv[puniq[(ifinal_l,l0)]]
                                                    if (abs(Fil[i0,l0])<smallb): continue
                                                    if (abs(Flk[l0,k0])<smallb): continue
                                                    contr = -Fij[i0,j0]*Fjk[j0,k0]*Flk[l0,k0]*Fil[i0,l0]
                                                    akey = (ju,ku,lu)
                                                    if akey in OCAs[iu]:
                                                        OCAs[iu][akey] += contr
                                                    else:
                                                        OCAs[iu][akey] = contr
        # OCA diagramms are renumbered to (i0,i1,i2,i3,ib1,ib2) where
        # i0,i1,i2,i3 are atomic states involved in the diagram and ib1,ib2 are the two bath
        # propagators in the diagram.
        OCAf={}
        for i in range(len(iequiv)):
            for ib1 in range(len(kbth)):
                for ib2 in range(len(kbth)):
                    for ks in list(OCAdiag[ib1][ib2][i].keys()):
                        if abs(OCAdiag[ib1][ib2][i][ks])<1e-10: continue
                        if (ib2<=ib1):
                            new_key = (i,) + ks + (ib1,ib2)
                            OCAf[new_key] = OCAdiag[ib1][ib2][i][ks]
                        else: # Due to time invariance, some diagrams are equivalent. For example
                              #  0 (b1) 1 (b2) 4 (b1) 2 (b2) 0    and      0 (b2) 2 (b1) 4 (b2) 1 (b1) 0
                            new_key = (i,) + (ks[2],ks[1],ks[0]) + (ib2,ib1)
                            OCAf[new_key] = OCAdiag[ib1][ib2][i][ks]

    # Here we regroup N_a =F^+ F  in more convenient way
    rFpF = zeros((len(iequiv),len(kbth)))
    for i in range(len(equiv)):
        df = FpF[0]-FpF[i]
        #for j in range(len(df)):  # if we don't do that, the occupancy can be annoyingly negative
        #    if abs(df[j])<1e-10: df[j]=0
        rFpF[equiv[i]] += df
    for i,ii in enumerate(iequiv):
        rFpF[i]/=len(iequiv[i]) # it has to be average, not the sum

    # Bubble contains diagrams named b (back). We need to construct from these also the other diagramms
    # which are called f (forward).
    forw_bubs = CreateEmpty2D_Dict(len(iequiv), len(kbth))
    for i in range(len(iequiv)):
        for b in range(len(kbth)):
            for ks in back_bubs[i][b]:
                forw_bubs[ks][b][i]= back_bubs[i][b][ks]

    # We want to have energy of each OCA-pseudoparticle ready
    Eq=zeros(len(iequiv))
    Egs=zeros(len(iequiv))
    Nq=zeros(len(iequiv),dtype=int)
    Nc0=0; Eg0=ipE[0]
    for i,ii in enumerate(iequiv):
        Nq[i] = ipN[ii[0]]
        Eq[i] = ipE[ii[0]]
        Egs[i] = Eg0   # ground state in this sector
        if (Nq[i]>Nc0):
            Nc0=Nq[i]
            Eg0=Eq[i]
            Egs[i]=Eg0

    # last renumbering for printing!
    # Only occupancies in n=[....] need to be kept. The rest of the diagrams is ignored.
    pu=-ones(len(iequiv),dtype=int)
    pl=0
    for i,ii in enumerate(iequiv):
        if Nq[i] in n:
            pu[i]=pl
            pl+=1

    # Printing output for OCA
    foca=open('cix.dat', 'w')
    print('# Input file for OCA impurity solver.', 'nlevels=', nlevels, 'Eimp=', Eimp, 'mOCA=', mOCA,  'Eoca=', Eoca, file=foca)
    print(len(kbth), (("%d "*len(kbth)) % tuple(map(len,kbth))), pl, 0, end=' ', file=foca)
    print('# Number of impurity levels, their degeneracy and number of local valence and local core states', file=foca)
    print('%3s' % '#', end=' ', file=foca)
    print(("%6s")*len(kbth) % tuple(['N'+str(x) for x in range(len(kbth))]), end=' ', file=foca)
    print("%4s" % 'Mtot', '%4s' % 'deg', '%12s' % 'Eatom', end=' ', file=foca)
    print(("%3s"%'#b')*len(back_bubs[i]), end=' ', file=foca)
    print(("%3s"%'#f')*len(forw_bubs[i]), file=foca)

    # write total spin S of PPs to file
    #ppfile=open('ppspin.dat', 'w')    # no need ppspin.dat, by Mateus

    for i,ii in enumerate(iequiv):
        if pu[i] <0 : continue  # state not used

        #print("%3d" % pu[i], "   ", ipS[ii[0]], file=ppfile)    # no need ppspin.dat, by Mateus

        print("%3d" % pu[i], (("%6.2f")*len(kbth) % tuple(rFpF[i])), "%4d" % Nq[i], end=' ', file=foca) #ipN[ii[0]],
        print("%4d" % len(ii), end=' ', file=foca)

        Eatom = ipE[ii[0]]  # This is the atomic energy
        # The single particle energy part will be added back inside the impurity solver,
        # therefore the energy should be subtracted
        for ib in range(len(kbth)): Eatom -= rFpF[i][ib]*Ed[ib]

        print("%12.6f" % Eatom, end=' ', file=foca)

        for b in range(len(kbth)): # delete diagrams which include states that were removed
            for ks in list(back_bubs[i][b].keys()):
                if pu[ks]<0: # this diagram involves state not considered
                    del back_bubs[i][b][ks]
            for ks in list(forw_bubs[i][b].keys()):
                if pu[ks]<0: # this diagram involves state not considered
                    del forw_bubs[i][b][ks]

        print(("%3d"*len(back_bubs[i])) % tuple(map(len,back_bubs[i])), end=' ', file=foca)
        print(("%3d"*len(forw_bubs[i])) % tuple(map(len,forw_bubs[i])), end=' ', file=foca)
        print('  ', end=' ', file=foca)

        for b in range(len(kbth)):
            for ks in back_bubs[i][b]:
                print("%10.6f x %-3d  " % (back_bubs[i][b][ks], pu[ks]), end=' ', file=foca)

        for b in range(len(kbth)):
            for ks in forw_bubs[i][b]:
                print("%10.6f x %-3d  " % (forw_bubs[i][b][ks], pu[ks]), end=' ', file=foca)

        print('# S ', ipS[ii[0]], file=foca) #, ' Eatom=', ipE[ii[0]]

    if qOCA:
        print('# OCA diagrams, information is (pp0,pp1,pp2,pp3) (b1,b2) fact , where pp is pseudoparticle and b is bath', file=foca)
        OCAF = list(OCAf.items())
        OCAF.sort(key=cmp_to_key(lambda x,y: cmp(y[1],x[1])))   # by Mateus
        for i in range(len(OCAF)):
            excitedE = [Eq[j]-Egs[j] for j in OCAF[i][0]]
            states_involved = [pu[l] for l in OCAF[i][0][:4]]
            # One of the states is not considered
            if (-1 in states_involved): continue
            # We take it into account only if all states that are involved,
            # have energy close to the ground state energy for this occupancy
            if max(excitedE)>Eoca:  continue
            # Matrix element negligible
            if abs(OCAF[i][1])<mOCA: continue
            if not (Nq[OCAF[i][0][1]] in Ncentral): continue

            print("%3d %3d %3d %3d   " % tuple(states_involved), end=' ', file=foca) #tuple([pu[l] for l in OCAF[i][0][:4]]),
            print("%2d %2d" % tuple(OCAF[i][0][4:]), end=' ', file=foca)
            print(OCAF[i][1], end=' ', file=foca)
            print('   #', [Eq[j]-Egs[j] for j in OCAF[i][0]], file=foca)

    ############################################################################################
    ######## End of the part which generates input for OCA solver                        #######
    ############################################################################################



