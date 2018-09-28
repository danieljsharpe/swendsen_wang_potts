import numpy as np
import random
import hoshen_kopelman_2

N = 4 # length of array

latt = np.zeros((N,N))

for i in range(0,N):
    for j in range(0,N):
        latt[i,j] = random.randint(0,1)

#latt = np.array([[1,0,1,1],[0,0,1,0],[0,1,1,1],[1,0,0,0]])
#latt = np.array(([0,1,0,1],[0,1,1,1],[0,0,1,0],[1,1,1,1]))
#latt = np.array(([1,0,1,1],[0,0,1,0],[0,1,0,1],[0,1,0,0]))


#latt = np.array(([1,0,1,0],[1,1,1,0],[0,1,0,1],[1,0,1,0]))
#latt = np.array(([0,1,1,0],[1,1,0,0],[0,0,0,1],[1,1,1,1]))


clust = hoshen_kopelman_2.hoshkop(latt)

print latt
print clust
