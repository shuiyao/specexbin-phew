import matplotlib.pyplot as plt
from scipy import array, log10, linspace, histogram
import ioformat

# f1 = "specztau.p50n288o5.37.200_0.vpm"
# f2 = "specztauw.p50n288o5pp.37.200_0.vpm"
f1 = "./autovp/HI.vpm"
f2 = "./autovp/OVI.vpm"

NHI1, v1 = ioformat.rcol(f1, [1,2])
for i in range(len(NHI1)): NHI1[i] = log10(NHI1[i]) + 13
NHI2, v2 = ioformat.rcol(f2, [1,2])
for i in range(len(NHI2)): NHI2[i] = log10(NHI2[i]) + 13
plt.plot(v1, NHI1, "b.")
plt.plot(v2, NHI2, "r.")
# x, y = histogram(NHI1, bins=linspace(13., 19., 20))
# plt.plot(x, y[1:], "b-")
# x, y = histogram(NHI2, bins=linspace(13., 19., 20))
# plt.plot(x, y[1:], "r-")
plt.show()
