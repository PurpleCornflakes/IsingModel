import matplotlib.pyplot as plt
import numpy as np

# mdata = np.loadtxt("m_T", delimiter = ",")
edata = np.loadtxt("e_T", delimiter = ",")

# m_ave = mdata[:,0]
# mm = mdata[:,1] # susceptiblity
e_ave = edata[:,0]
c = edata[:,1] # capacity

T = np.arange(2.0, 3.5, 0.2)

# plt.plot(T, e_ave/32/32)
# plt.title("M_T")
# plt.title("chi_T")
# plt.title("e_T")
plt.title("ee_T")

# plt.plot(T, e_ave)
plt.plot(T, c*T**2)
# plt.savefig("e_T.png", dpi = 200)
plt.savefig("ee_T.png", dpi = 200)
plt.show()
