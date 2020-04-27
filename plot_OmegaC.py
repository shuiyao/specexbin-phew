import matplotlib.pyplot as plt
from scipy import linspace

tauname = "specaim.p50n288ezw15.4_10.98_8.680.x+25kpc"
binname = "binzfile.p50n288ezw15.4_10.98_8.680.x+25kpc"

fbin = open(binname, "r")
zlist, mlist, Clist = [], [], []
for line in fbin:
    if (line[0] != "#"):
        spt = line.split()
        mlist.append(float(spt[2]))
        zlist.append(float(spt[1]))
        Clist.append(float(spt[6]))
fbin.close()

zbins = linspace(0.26, 0.24, 30.)
i, j = 0, 0
Cmass = 0.
totmass = 0.
Omega_C = []
for z in zlist:
    if z < 0: break
    if z < zbins[i]:
        if totmass == 0.:
            Omega_C.append(0.)
        else:
            Omega_C.append(Cmass/totmass)
        Cmass = 0.
        totmass = 0.
        i = i + 1
    else:
        totmass += mlist[j]
        Cmass += mlist[j] * Clist[j]
    j = j + 1

plt.plot(zbins[1:], Omega_C, "o-")
plt.yscale("log")
plt.show()
