import matplotlib.pyplot as plt

#f = open("tautest2.txt", "r")
#f = open("specaim.p50n288ezw15.4_10.98_8.680.x+25kpc", "r")
#f = open("specaim.p50n288ezw15.70_10.77_8.368.x+25kpc", "r")
# f = open("specaimc.test.halo2", "r")
f = open("specztauw.p50n288o5pp.37.200_0", "r")
# f = open("ztau.37.5_0", "r")
mode = 0 # Optical depth
#mode = 1 # rho, T

#f = open("specaim.test.halo", "r")
#f = open("specaim.test.halo1.original", "r")

# cols = f.readline().split()
# print len(cols)
# taus = []
# i = 7
# while(i<43):
#     taus.append(cols[i])
#     i = i + 4
# print taus

z, tauHI, tauOVI = [], [], []
rho, T, MZ = [], [], []
for line in f:
    z.append(float(line.split()[0]))
    rho.append(float(line.split()[1]))
    T.append(float(line.split()[2]))
    MZ.append(float(line.split()[3]))
    tauHI.append(float(line.split()[7]))
    tauOVI.append(float(line.split()[19]))
if(mode == 0):
    plt.plot(z, tauHI)
    plt.plot(z, tauOVI, "--")
    plt.yscale("log")
if(mode == 1):
    plt.plot(z, T, "--")
    plt.plot(z, rho)

#plt.axis([0.24, 0.26., 0.1, 1.e8])
f.close()
