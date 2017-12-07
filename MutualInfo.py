import numpy as np 
import matplotlib.pyplot as plt 
from matplotlib import pylab
from sklearn.metrics.cluster import normalized_mutual_info_score, adjusted_mutual_info_score


ising_slice = np.loadtxt("slice_Tc", delimiter = ",")
adpNN_r = np.load("r1_64.npy")

def plot_MI_distance(data, ax, marker, dmax = 25, method = 'NMI'):
	if len(data.shape) == 1:
		assert(len(data)%dmax == 0)
		data = data.reshape(len(data)//dmax,dmax)
	N,dmax = data.shape
	print(dmax)
	if method == 'NMI':
		NMIs = np.zeros(dmax)
		for d in range(1,dmax):
			NMIs[d] = normalized_mutual_info_score(data[:,0], data[:,d])
		MIs = NMIs

	if method == 'AMI':
		AMIs = np.zeros(dmax)
		for d in range(1,dmax):
			AMIs[d] = adjusted_mutual_info_score(data[:,0], data[:,d])
		MIs = AMIs

	if method == 'self_MI':
		self_MIs = np.zeros(dmax)
		for d in range(1,dmax):
			self_MIs[d] = mutual_info(data[:,0], data[:,d])
		MIs = self_MIs
	# return ax.plot(range(1,dmax), MIs[1:dmax],'o')	
	return ax.plot(np.log10(range(1,dmax)), np.log10(MIs[1:dmax]),marker)[0]

def plot_all():
	plt.close('all')
	axes = [0]
	fig,axes[0] = plt.subplots(1)
	fig.set_size_inches(10, 5)

	ising_NMIs = plot_MI_distance(data=ising_slice[:,:25], ax=axes[0], marker='ro')
	ising_AMIs = plot_MI_distance(data=ising_slice[:,:25], ax=axes[0], method='AMI', marker='r-')
	ising_MIs = plot_MI_distance(data=ising_slice[:,:25], ax=axes[0], method='self_MI', marker='r--')
	adpNN_NMIs = plot_MI_distance(data=adpNN_r, ax=axes[0], marker='bo')
	adpNN_AMIs = plot_MI_distance(data=adpNN_r, ax=axes[0], method='AMI', marker='b-')
	adpNN_MIs = plot_MI_distance(data=adpNN_r, ax=axes[0], method='self_MI', marker='b--')

	plt.legend((ising_NMIs, ising_AMIs, ising_MIs, adpNN_NMIs, adpNN_AMIs, adpNN_MIs),\
		('Ising_NMI', 'Ising_AMI', 'Ising_MIs', 'AdpNN_r_NMI', 'AdpNN_r_AMI', 'AdpNN_MIs'))
	plt.show() 

def H(Px):
	return np.sum(-Px * np.log2(Px))

def mutual_info(x, y, bins = [[-1.5, 0.5, 1.5],[-1.5, 0.5, 1.5]]):
	Nxy = np.histogram2d(x, y, bins = bins)[0]
	Nx = np.sum(Nxy, axis=0)
	Ny = np.sum(Nxy, axis=1)
	Pxy = Nxy/np.sum(Nxy)
	Px = Nx/np.sum(Nx)
	Py = Ny/np.sum(Ny)
	mi = H(Px) + H(Py) - H(Pxy)
	return mi/np.sqrt(H(Px)*H(Py))
 
# mutual_info(ising_slice[:,0], ising_slice[:,1])

plot_all()