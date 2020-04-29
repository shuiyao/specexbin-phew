import matplotlib.pyplot as plt

ion="MgII"
fname = ion+".cln"
wvlen, v, flux = [], [], []
f = open(fname, "r")
for line in f:
    spt = line.split()
    wvlen.append(float(spt[0]))
    v.append(float(spt[1]))
    flux.append(float(spt[2]))
f.close()
plt.plot(v, flux)
