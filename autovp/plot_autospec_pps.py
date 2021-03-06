import matplotlib.pyplot as plt
from scipy import array
import matplotlib as mpl

ions=["HI", "MgII","CIV","OVI"]
clrs=["black", "blue", "green", "red"]
i = 0

pfs = []
voffset = 300.

i = 0
ps = []
for ion in ions:
    fname = "pp."+ion+".cln"
    wvlen, v, flux = [], [], []
    f = open(fname, "r")
    for line in f:
        spt = line.split()
        wvlen.append(float(spt[0]))
        v.append(float(spt[1])+voffset*i)
        flux.append(float(spt[2]))
    f.close()
    p, = plt.plot(v, flux, "-", color=clrs[i])
    pf, = plt.fill(v+v[0:1], flux+flux[0:1], color=clrs[i], alpha=0.2)
    if(ion == "MgII"): pfs.append(pf)
    ps.append(p)
    i = i + 1

i = 0
for ion in ions:
    fname = "ppe."+ion+".cln"
    wvlen, v, flux = [], [], []
    f = open(fname, "r")
    for line in f:
        spt = line.split()
        wvlen.append(float(spt[0]))
        v.append(float(spt[1])+voffset*i)
        flux.append(float(spt[2]))        
        # flux.append(1.0-float(spt[2]))
    f.close()
    plt.plot(v, flux, "-", color=clrs[i], linewidth=2)
    pf, = plt.fill(v+v[0:1], flux+flux[0:1], color="none", hatch="...", ec=clrs[i], alpha=0.5)
    # pf, = plt.fill(v+v[0:1], flux+flux[0:1], color=clrs[i], alpha=0.2)
    pf = plt.Rectangle((v[0], flux[0]), 0, 0, fill=False, color="blue", hatch="....")
    if(ion == "MgII"): pfs.append(pf)
    i = i + 1
    
# plt.axis([66600., 67200., 0.0, 1.0])
plt.axis([66600., 67200.+voffset*3, 0.0, 1.0])
l1 = plt.legend(ps, ions, loc=3)
l2 = plt.legend(pfs, ["Smoothed", "P-by-P"], loc=4, frameon=1)
frame = l2.get_frame()
frame.set_facecolor('none')
frame.set_edgecolor('black')
plt.gca().add_artist(l1)
plt.xlabel("v [km/s]")
plt.ylabel("Flux")
plt.show()
