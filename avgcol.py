from scipy import linspace, sqrt, pi
from scipy.integrate import quad
import matplotlib.pyplot as plt

def cs32(u):
    C1 = 2.546479089470
    C2 = 15.278874536822
    C5 = 5.092958178941
    if (u < 0.): return 0.
    elif (u < 0.5): return C1 + C2 * (u - 1.0) * u * u
    elif (u < 1.0): return C5 * (1.0 - u)**3
    else: return 0.

def intz_cs32(z, *args):
    r = sqrt(z * z + args[0] * args[0])
    if(r > 1.0): r = 1.0
    return cs32(r)
    
def int_cs32(u):
    return cs32(u) * 4. * pi * u * u

nbins = 1000
d = linspace(0.0, 1.0, nbins)
intr = []
halfdz = []
totint = 0.0
C = 10.
NCloud = 1000.0 # Each cloud has mass 1./Ncloud
nclouds, colclouds = [], []

for i in range(len(d)):
    zmin = -sqrt(1.0 - d[i] * d[i])
    zmax = -zmin
    halfdz.append(zmax)
    intr.append(quad(intz_cs32, zmin, zmax, args=(d[i],))[0])
    # rho * intr is the column density
    # Assuming rhoa = 1.0; rhoc = C * rhoa
    # integrate over all columns gives 1.0
    # define clumping factor as C = rhoc / rhoa > 1.0
    #    therefore Rc / h = C**(-1./3.)*(Ncloud)**(-1./3.) < 1.0
    nclouds.append(intr[-1]*C**(-2./3.)*NCloud**(1./3.))
    avg_dzcloud = 4./3.*C**(-1./3.)*NCloud**(-1./3.)
    colclouds.append(nclouds[-1]*avg_dzcloud*C)
    totint += 2.*d[i]*pi*1.0/float(nbins)*intr[-1]

plt.plot(d, intr, "b-")
plt.plot(d, halfdz, "g-")
plt.plot(d, colclouds, "r-")
plt.plot(d, nclouds, "r--")
print totint
plt.show()
