#Hoshen-Kopelman algorithm if the array arr does not have periodic boundary conditions applied

import numpy as np

def hoshkop(arr):
    nclust = 0
    labels = [i for i in range(len(arr)**2)]
    clust = np.zeros((len(arr),len(arr)))
    for i in range(0,len(arr)):
        for j in range(0,len(arr)):
            if arr[i,j] == 1:
                if i == 0 and j == 0: # top left cell
                    nclust += 1
                    clust[i,j] = nclust
                elif i == 0 and j != 0: # do not check above
                    if arr[i,j-1] == 0:
                        nclust += 1
                        clust[i,j] = nclust
                    else:
                        clust[i,j], labels = find(clust[i,j-1],labels)
                elif i != 0 and j == 0: #do not check to left
                    if arr[i-1,j] == 0:
                        nclust += 1
                        clust[i,j] = nclust
                    else:
                        clust[i,j], labels = find(clust[i-1,j],labels)
                else:
                    if arr[i-1,j] == 0 and arr[i,j-1] == 0:
                        nclust += 1
                        clust[i,j] = nclust
                    elif arr[i-1,j] != 0 and arr[i,j-1] == 0:
                        clust[i,j], labels = find(clust[i-1,j],labels)
                    elif arr[i-1,j] == 0 and arr[i,j-1] != 0:
                        clust[i,j], labels = find(clust[i,j-1],labels)
                    else:
                        labels = union(clust[i,j-1],clust[i-1,j],labels)
                        clust[i,j], labels = find(clust[i-1,j],labels)
    clust = relabel(clust,labels)
    return clust

def find(m,labels):
    y = int(m)
    while (labels[y] != y):
        y = labels[y]
    while (labels[int(m)] != int(m)):
        z = labels[int(m)]
        labels[int(m)] = y
        m = z
    return y, labels

def union(m,n,labels):
    labels[int(find(m,labels)[0])] = int(find(n,labels)[0])
    return labels

def uf_makeset():
    return

def relabel(clust,labels):
    newlabels = np.zeros(len(clust)*2)
    for i in range(0,len(clust)):
        for j in range(0,len(clust)):
            if clust[i,j] != 0:
                x, labels = find(clust[i,j],labels)
                if newlabels[x] == 0:
                    newlabels[0] += 1
                    newlabels[x] = newlabels[0]
                clust[i,j] = newlabels[x]
    return clust
