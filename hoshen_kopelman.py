#module containing functions to apply Hoshen-Kopelman algorithm (i.e. union-find type algorithm) to an array
#for use with swendsenwang_potts, i.e. is for the mapping of a site percolation problem onto a bond percolation problem
#The algorithm accounts for periodic boundary conditions

import numpy as np

#scan and labelling clusters numbers on array 'latt'
def hoshkop(latt,bonds):
    nclust = 0 # no. of clusters found (not accounting for unions)
    labels = [i for i in range(len(latt)**2)]
    clust = np.zeros((len(latt),len(latt)), dtype=int)
    # first cycle - ignore PBCs
    for i in range(0,len(latt)):
        for j in range(0,len(latt)):
            if i == 0 and j == 0: # first cell
                nclust += 1
                clust[i,j] = nclust
            elif i == 0 and j != 0: # first row - do not check vertical bonds
                if bonds[i,2*j] == 1: # horizontal bond i,j -> i,j-1
                    clust[i,j], labels = find(clust[i,j-1], labels)
                else:
                    nclust += 1
                    clust[i,j] = nclust
            elif i != 0 and j == 0: # first column - do not check horizontal bonds
                if bonds[i,2*j+1] == 1: # vertical bond i,j -> i-1,j
                    clust[i,j], labels = find(clust[i-1,j], labels)
                else:
                    nclust += 1
                    clust[i,j] = nclust
            else:
                if bonds[i,2*j] == 1 and bonds[i,2*j+1] == 0: # there exists a 'horizontal' bond i,j -> i,j-1 only
                    clust[i,j], labels = find(clust[i,j-1], labels)
                elif bonds[i,2*j] == 0 and bonds[i,2*j+1] == 1: # there exists a 'vertical' bond, i,j -> i-1,j only
                    clust[i,j], labels = find(clust[i-1,j], labels)
                elif bonds[i,2*j] == 1 and bonds[i,2*j+1] == 1: # there exists both a vertical and horizontal bond   
                    labels = union(clust[i,j-1], clust[i-1,j], labels)
                    clust[i,j], labels = find(clust[i-1,j], labels)
                else: # no bonds
                    nclust += 1
                    clust[i,j] = nclust
    clust, labels = check_pbc(clust, labels, bonds)
    cluster_sets = relabel(clust, labels)
    return cluster_sets

def union(m,n,labels):
    labels[int(find(m,labels)[0])] = int(find(n,labels)[0])
    return labels

# x = find(x) has the property: labels[x] == x, which is the defining property for x to be the representative member of its
# equivalence class
def find(m,labels):
    y = int(m)
    while labels[y] != y:
        y = labels[y]
    while (labels[int(m)] != int(m)):
        z = labels[int(m)]
        labels[int(m)] = y
        m = z
    return y, labels

def check_pbc(clust, labels, bonds):
    for j in range(np.shape(clust)[0]): # check top row
        if bonds[0,2*j] == 1 and bonds[0,2*j+1] == 0: # there exists a 'horizontal' bond i,j -> i,j-1 only
            clust[0,j], labels = find(clust[0,j-1], labels)
        elif bonds[0,2*j] == 0 and bonds[0,2*j+1] == 1: # there exists a 'vertical' bond, i,j -> i-1,j only
            clust[0,j], labels = find(clust[-1,j], labels)
        elif bonds[0,2*j] == 1 and bonds[0,2*j+1] == 1: # there exists both a vertical and horizontal bond   
            labels = union(clust[0,j-1], clust[-1,j], labels)
            clust[0,j], labels = find(clust[-1,j], labels)
        else: # no bonds
            continue
    for i in range(1, np.shape(clust)[0]): # check leftmost column (except top left cell)
        if bonds[i,0] == 1 and bonds[i,1] == 0: # there exists a 'horizontal' bond i,j -> i,j-1 only
            clust[i,0], labels = find(clust[i,-1], labels)
        elif bonds[i,0] == 0 and bonds[i,1] == 1: # there exists a 'vertical' bond, i,j -> i-1,j only
            clust[i,0], labels = find(clust[i-1,0], labels)
        elif bonds[i,0] == 1 and bonds[i,1] == 1: # there exists both a vertical and horizontal bond   
            labels = union(clust[i,-1], clust[i-1,0], labels)
            clust[i,0], labels = find(clust[i-1,0], labels)
        else: # no bonds
            continue
    return clust, labels

def relabel(clust, labels):
    cluster_sets = {}
    for i in range(np.shape(clust)[0]):
        for j in range(np.shape(clust)[0]):
            clust[i,j] = find(clust[i,j], labels)[0]
            while True:
                try:
                    cluster_sets[int(clust[i,j])].append((i,j))
                except KeyError:
                    cluster_sets[int(clust[i,j])] = []
                    continue
                break
    return cluster_sets
