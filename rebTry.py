import rebound
import numpy as np
import astropy.units as u
import matplotlib.pyplot as plt

kep33 = rebound.Simulation() # This sets up our system
# kep33.integrator = 'whfast'

# hmmm what units?

kep33.add(m=1) #star
kep33.add(m=1e-6, P=21.0*60*60*24) 
kep33.add(m=1e-6, P=30.0*60*60*24)
kep33.move_to_com() # Moving our system to the centre of mass frame, which should be in the star.


N = 30
times = np.sort(50*60*60*24*np.random.random(N)) # 30 randomly spaced observations
p1, p2 = np.zeros((N,2)), np.zeros((N,2))

for i, t in enumerate(times):
    kep33.integrate(times[i])

    p1[i] = np.array([kep33.particles[1].x, kep33.particles[1].y])
    p2[i] = np.array([kep33.particles[2].x, kep33.particles[2].y])

plt.scatter(p1[:, 0], p1[:, 1])
plt.scatter(p2[:, 0], p2[:, 1])
plt.show()



#print(kep33.status())

# plt.plot(kep33)
# plt.show()
