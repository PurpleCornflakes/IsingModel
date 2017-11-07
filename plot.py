import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("m_T", delimiter = ",")

m_ave = data[:,0]
mm = data[:,1] # susceptiblity
T = np.arange(2.0, 3.5, 0.2)

# plt.plot(T, e_ave/32/32)
plt.title("M_T")
# plt.title("chi_T")

plt.plot(T, m_ave)
# plt.plot(T, mm)
plt.savefig("m_T.png", dpi = 200)
plt.show()
# plt.savefig("chi_T.png")