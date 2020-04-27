import matplotlib.pyplot as plt
from scipy import log10, array

#f = open("tautest2.txt", "r")
#f = open("specaim.p50n288ezw15.4_10.98_8.680.x+25kpc", "r")
#f = open("specaim.p50n288ezw15.70_10.77_8.368.x+25kpc", "r")
# f = open("specaimc.test.halo2", "r")
#f = open("specztau.p50n288o5pp.37.5_0", "r")
mode = 0 # Optical depth
mode = 2 # Optical depth HI + T
#mode = 1 # rho, T
# fname1 = "ztau.37.5_0"
# fname2 = "ztauw.37.5_0"
#fname1 = "specztau.p50n288o5sm.37.2_0"
#fname2 = "ztauw.37.2_0"
# fname1 = "specztau.p50n288o5.10.200_0"
# fname2 = "specztauw.p50n288o5pp.10.200_0"

fname1 = "shortlos.p50n288o5.10kpc"
fname2 = "shortlos.p50n288o5pp.10kpc"

# fname1 = "p25.ztau.37.2_0"
# fname2 = "specztauw.p25n144gwl.37.2_0"

f = open(fname1, "r")
z, tauHI, tauOVI = [], [], []
rho, T, MZ, logT = [], [], [], []
for line in f:
    z.append(float(line.split()[0]))
    rho.append(float(line.split()[1]))
    logT.append(float(line.split()[2]))
    MZ.append(float(line.split()[3]))
    tauHI.append(log10(float(line.split()[7])))
    tauOVI.append(float(line.split()[19]))
if(mode == 0):
    plt.plot(z, tauHI, "b-")
    plt.plot(z, tauOVI, "b:")
    plt.yscale("log")
if(mode == 1):
    plt.plot(z, T, "b:")
    plt.plot(z, rho, "b-")
if(mode == 2):
    plt.plot(z, array(logT)-3.0, "b:")
    plt.plot(z, rho, "b--")            
    plt.plot(z, tauHI, "b-")
f.close()
print sum(tauHI)

f = open(fname2, "r")
z, tauHI, tauOVI = [], [], []
rho, T, MZ, logT = [], [], [], []
for line in f:
    z.append(float(line.split()[0]))
    rho.append(float(line.split()[1]))
    logT.append(float(line.split()[2]))
    MZ.append(float(line.split()[3]))
    tauHI.append(log10(float(line.split()[7])))
    tauOVI.append(float(line.split()[19]))
if(mode == 0):
    plt.plot(z, tauHI, "r-")
    plt.plot(z, tauOVI, "r:")
    plt.yscale("log")
if(mode == 1):
    plt.plot(z, T, "r:")
    plt.plot(z, rho, "r-")
if(mode == 2):
    plt.plot(z[1:], array(logT[:-1])-3.0, "r:")
    plt.plot(z[1:], rho[:-1], "r--")        
    plt.plot(z[1:], tauHI[:-1], "r-")

f.close()    
print sum(tauHI)

plt.show()    
f.close()
