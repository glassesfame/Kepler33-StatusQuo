import os
import rebound
import numpy as np
import pandas as pd
from scipy import stats
import astropy.constants as ast
import matplotlib.pyplot as plt


def fCalc(b4T0, P, e):
    '''
    To approximate the starting positions. First, we find the mean ananomaly, taking T0 as periapsis.
    Then, we approximate the eccentric ananomaly as the mean ananomaly as e is small. We then directly
    derive our true ananomaly from our eccentric ananomaly, given as the alternate form from
    https://en.wikipedia.org/wiki/True_anomaly
    '''
    M = (P - b4T0)/P # mean ananomaly approximately eccentric ananomaly when e is small.
    # The calculation to derive the true 
    beta = e/(1+np.sqrt(1-e**2))
    temp = beta*np.sin(M)/(1-beta*np.cos(M))
    f = M + 2*np.arctan(temp)
    return f

    
# Handling our data which are the parameters for the simulation. 
exoNames = ['b', 'c', 'd', 'e', 'f']
df = pd.read_csv('kep33Data(Sheet1).csv', names=exoNames)

# Getting the required data for inputting into the parameters.
getData = lambda idx: np.array(df.iloc[idx], dtype=float)
eArr = getData(1)
fCol = fCalc(getData(4), getData(0), eArr) # true anomaly gives us our initial position but this must be calculated from our timings.

kep33 = rebound.Simulation() # This sets up our system
# kep33.integrator = 'whfast'
# Our star is at the centre of the system with no velocity etc, so no orbital elements for it except mass.
kep33.add(m=1.29) 
# Populating the simulation with our planets right now. 
for i, colName in enumerate(exoNames):
    kep33.add(m=(getData(3)* (ast.M_earth.value/ast.M_sun.value))[i], P=getData(0)[i], e=getData(1)[i], inc=np.pi/2, f=fCol[i])
    
kep33.move_to_com() # Moving our system to the centre of mass frame, which should be in the star.

times = np.arange(0, 1600, 2)
N = len(times) 
pTot = np.zeros((N,6,2))

for i, t in enumerate(times):
    kep33.integrate(t)
    for j in range(pTot.shape[1]):
        pTot[i, j] = np.array([kep33.particles[j].x, kep33.particles[j].y])

fig, ax = plt.subplots(figsize=(10,8))
for i, exo in enumerate(exoNames):
    ax.scatter(pTot[:, i+1, 0], pTot[:, i+1, 1], label=exo, s=1)
ax.legend()
plt.show()