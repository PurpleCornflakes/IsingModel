import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("m_T", delimiter = ",")

m_ave = data[:,0]
mm = data[:,1] # susceptiblity
T = np.arange(2.0, 3.5, 0.2)

# plt.plot(T, e_ave/32/32)
plt.title("m distribution")
# plt.title("chi distribution")

plt.plot(T, m_ave)
# plt.plot(T, mm)
plt.show()
plt.savefig("m_T.png")
# plt.savefig("chi_T.png")