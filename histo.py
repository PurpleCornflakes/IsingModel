import matplotlib.pyplot as plt
import numpy as np

# ls = plt.imread("ising100T15.png")
# cs = plt.imread("ising100Tc.png")
# hs = plt.imread("ising100T5.png")
# imgs = [ls, cs, hs]
# titles = ["low T", "critical", "high T"]
# fig0, axes = plt.subplots(ncols = len(imgs), figsize = [14,6])

# for i, ax in enumerate(axes):
#     ax.imshow(imgs[i])
#     ax.set_title(titles[i])
#     ax.set_xticks([])
#     ax.set_yticks([])

# 100*100
# iteration 5000
l = np.loadtxt("m_dist_T15")
c = np.loadtxt("m_dist_Tc")
h = np.loadtxt("m_dist_T50")

# bins = np.arange(-5000,5000,500)


fig, axes = plt.subplots(ncols=3)
axes[0].hist(l,color = "blue")
axes[1].hist(c, color = "cyan")
axes[2].hist(h, color = "red")


# fig.set_title("100*100 spins @kT = 1.5, 2.26 and 5.0")
# fig.set_xlabel("M")
# fig.set_ylabel("counts")
plt.title("M_distribution")
plt.tight_layout()
plt.savefig("M_distribution.png", dpi = 200)
plt.show()