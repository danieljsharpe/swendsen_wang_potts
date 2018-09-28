# Script to simulate q-state Potts glass using the Swendsen-Wang algorithm
# Daniel J. Sharpe
# 25.11.2017

import numpy as np
import matplotlib.pyplot as plt
import random
import pbc
import hoshen_kopelman
import analysis
from matplotlib import colors as mcolors
from matplotlib import colorbar
#from corrfunc import *

class Sw_qpotts(object):

    def __init__(self,J,T,q,N,nsteps,seed=17):
        self.J = J # interaction energy
        self.T = T # reduced temperature
        self.q = q # no. of possible spin states s = 0,1,2..
        self.N = N # dimension of square lattice
        self.nsteps = nsteps # no. of simulation steps
        random.seed(seed)

    def statement(self):
        return "Running simulation with parameters:\nJ = %f, T = %f, q = %g, N = %g, nsteps = %g" \
              % (self.J, self.T, self.q, self.N, self.nsteps)

    def getlattice(self):
        latt = np.zeros((self.N,self.N))
        latt = pbc.Periodic_Lattice(latt) # apply PBCs
        for i in range(0,self.N):
            for j in range(0,self.N):
                latt[i,j] = random.randint(0,self.q-1)
        return latt

    def getbonds(self,bonds,latt):
        nbonds = 0
        for i in range(len(latt)):
            for j in range(0,len(latt)):
                if latt[i,j] == latt[i,j-1]:
                    bonds[i,2*j] = 1 # 'horizontal' bonds
                    nbonds += 1
                else:
                    bonds[i,2*j] = 0
                if latt[i,j] == latt[i-1,j]:
                    bonds[i,2*j+1] = 1 # 'vertical' bonds
                    nbonds += 1
                else:
                    bonds[i,2*j+1] = 0
        return bonds

    def del_bonds(self, bonds, latt):
        for i in range(len(bonds)):
            for j in range(2*len(bonds)):
                if bonds[i,j] == 1:
                    if self.J <= 0.0:
                        if random.uniform(0.0,1.0) <= np.exp(4.0*(self.J/self.T)):
                            bonds[i,j] = 0 # deleted bond
                    elif self.J > 0.0:
                        if random.uniform(0.0,1.0) <= np.exp(-4.0*(self.J/self.T)):
                            bonds[i,j] = 0 # deleted bond
        return bonds, latt


    def spinflip(self, latt, cluster_sets):
        for key in cluster_sets:
            newq = random.randint(0,self.q-1)
            for cell in cluster_sets[key]:
                latt[cell] = newq
        return latt

    def simulate_sw(self):
        latt = self.getlattice()
        bonds = np.zeros((self.N,2*self.N))
        print self.statement()
        for i in range(0,self.nsteps):
            bonds = self.getbonds(bonds,latt)
            bonds, latt = self.del_bonds(bonds, latt)
            cluster_sets = hoshen_kopelman.hoshkop(latt,bonds)
            latt = self.spinflip(latt, cluster_sets)
        return latt

    def plot(self,latt):
        fig, ax = plt.subplots(1,1, figsize =(8,6)) # create grid and color bar as separate 'plots'
#        plt.matshow(latt)
        cmap = plt.cm.jet # define the colormap
        cmaplist = [cmap(i) for i in range(cmap.N)] # extract all colours from the .jet map
#        cmaplist[0] = (.5,.5,.5,1.0) # choose the first colour in the list, can do others likewise
        cmap = cmap.from_list('Custom cmap', cmaplist, cmap.N) # create new map
        bounds = np.linspace(0,self.q,self.q+1) # define bins
        norm = mcolors.BoundaryNorm(bounds, cmap.N) # normalise
        grid = ax.matshow(latt)
        ax2 = fig.add_axes([0.92, 0.1, 0.03, 0.8]) # create and set axis for colorbar
        cb = colorbar.ColorbarBase(ax2, cmap=cmap, norm=norm, spacing='proportional', \
             ticks=bounds, boundaries=bounds, format='%1i', orientation='vertical')
        ax.set_title("%i-state Potts model snapshot" % int(np.amax(latt)+1))
        ax.axis('off')
        ax2.set_ylabel("Spin state q = 0,1,2..", size=12)
        plt.show()
        return

sw1 = Sw_qpotts(-0.25,1.0,3,128,100)
latt1 = sw1.simulate_sw()
sw1.plot(latt1)
