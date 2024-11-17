import os
import rebound
import numpy as np
import pandas as pd
from scipy import stats
import astropy.constants as ast
import matplotlib.pyplot as plt

const ={
    'G' :  6.67430e10-11, #SI
    'R_b' : 0.0677, #AU
    'R_c' : 0.1189 , #AU
    'R_d' : 0.1662, #AU
    'R_e' : 0.2138, #AU
    'R_f' : 0.2535, #AU
    'P_b' : 5.7, # days
    'P_c' : 13.2, # days
    'P_d' : 21.8, # days
    'P_e' : 31.8, # days
    'P_f' : 41,# days
    'M_b' : 3.68, #Earth
    'M_c' : 0.39, #Earth
    'M_d': 3.91, #Earth
    'M_e' : 5.57, #Earth
    'M_f': 9.6, #Earth
    'AU' : 1.496e11, #m
    'M_star' : 1.29, # solar mass
    'M_sun' : 1.988e30, # kg
    'day_to_sec' : 86400,
    'M_earth' : 5972e24 #kg
}

def timingData(path):
    '''
    Extracting the timing data.
    '''
    calArr = [] # the calculated timings based on the periods.
    t0Arr = [] # the actual timings measured by Kepler

    for file in os.listdir(path):
        # Making sure we only take the data files.
        if 'koi0707.' in file:
            d = np.loadtxt(file)
            
            calArr.append(d[:, 0])
            t0Arr.append(d[:, 0] + d[:, 1])
    
    # Find the period of the data, which should be the mode of the differences. 
    P = np.array([stats.mode(arr[1:] - arr[:-1], keepdims=False).mode for arr in calArr])
    # We sort by the period because the data files are not in order. 
    sortIdx = np.argsort(P)
        
    return np.array(calArr, dtype=object)[sortIdx], P[sortIdx]

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

def getF(calArr, P, e):
    '''
    Extracting the data necessary to calculate the true ananomaly which allows us to calculate
    the starting position. 
    '''
    # The time stamp which we will start our calculations at. 
    initT = [arr[0] for arr in calArr]
    firstT0 = initT[np.argmin(initT)]
    b4T0 = initT - firstT0
    return fCalc(b4T0, P, e)
    
# Handling our data which are the parameters for the simulation. 
df = pd.read_csv('kep33Data(Sheet2).csv')
exoNames = df.columns[1:-1]
predArr, period = timingData(None)
# Getting the required data for inputtinf into the parameters.
eArr = np.array(df.iloc[1][1:-1], dtype=float)
fCol = getF(predArr, period, eArr) # true anomaly gives us our initial position but this must be calculated from our timings.

kep33 = rebound.Simulation() # This sets up our system
# kep33.integrator = 'whfast'
# Our star is at the centre of the system with no velocity etc, so no orbital elements for it except mass.
kep33.add(m=1.29) 
# Populating the simulation with our planets right now. 
for i, colName in enumerate(exoNames):
    exoPars = np.array(df[colName],dtype=float)
    kep33.add(m=exoPars[-1]*ast.M_earth.value/ast.M_sun.value, P=period[i], e=exoPars[1], inc=np.pi/2, f=fCol[i])
    
kep33.move_to_com() # Moving our system to the centre of mass frame, which should be in the star.

times = np.arange(0, 1600, 2)
N = len(times) 
pTot = np.zeros((N,5,2))

for i, t in enumerate(times):
    kep33.integrate(t)
    for j in range(pTot.shape[1]):
        pTot[i, j] = np.array([kep33.particles[j+1].x, kep33.particles[j+1].y])

fig, ax = plt.subplots(figsize=(10,8))
ax.scatter(kep33.particles[0].x, kep33.particles[0].y)
print(kep33.particles[0].x, kep33.particles[0].y)
for i, exo in enumerate(exoNames):
    ax.scatter(pTot[:, i, 0], pTot[:, i, 1], label=exo, s=1)
ax.legend()
plt.show()