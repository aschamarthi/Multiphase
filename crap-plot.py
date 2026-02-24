# from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
# from mpl_toolkits.mplot3d import Axes3D
# import scipy.interpolate as interpolate


# # For latex font, i guess so 
plt.rc('text', usetex=True)
plt.rc('font', family='arial')
#Set global matplotlib parameters in script or in /home/$USER/.matplotlib/matplotlibrc
# plt.rcParams['axes.linewidth'] = 1.5
# plt.rcParams['xtick.major.size'] = 8
# plt.rcParams['xtick.minor.size'] = 4
# plt.rcParams['ytick.major.size'] = 6
# plt.rcParams['ytick.minor.size'] = 3
plt.rcParams.update({'font.size': 10})



# x,y,d1,d2,u,v,p,a,grad=np.loadtxt('Rslt0002.plt', delimiter=None, unpack=True,skiprows=3)
x,y,d1,u,v,p,a,grad,phi=np.loadtxt('Rslt0026.plt', delimiter=None, unpack=True,skiprows=3)

# x,y,phi,grad=np.loadtxt('new1.plt', delimiter=None, unpack=True,skiprows=2)

g=512*4
k=768*4


x = x.reshape(g,k)
y = y.reshape(g,k)
# d1 = d1.reshape(g,k)
# p = p.reshape(g,k)
# u = u.reshape(g,k)
# v = v.reshape(g,k)
a = a.reshape(g,k)
# phi = phi.reshape(g,k)
grad = grad.reshape(g,k)

# blah=np.log(abs(grad)+1)

# flop = np.exp(-.09*abs(grad))#/grad.max())

# print(blah.max())
# plt.imshow(phi, vmin = 0, vmax = 1.5, cmap=plt.cm.gray_r, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()])
# plt.contour(x,y,a,38,linewidths=0.1,colors=('k'))
# plt.contour(x,y,np.log(abs(grad)+1),40,cmap='jet')
plt.imshow(grad, vmin = 1, vmax = 2400, cmap=plt.cm.Blues, origin='upper', 
           extent=[x.min(), x.max(), y.min(), y.max()])
# plt.imshow(grad, vmin = 1, vmax = 75, cmap=plt.cm.gray_r, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()])
# plt.imshow(phi, vmin = 1, vmax = 1.5, cmap=plt.cm.gray_r, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()])
# plt.imshow(blah, vmin = 0, vmax = 3, cmap=plt.cm.gray_r, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()])
# plt.imshow(flop, vmin = 0.0, vmax = 1.0, cmap=plt.cm.gray, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()])
# fig, ax = plt.subplots()
# CS = ax.contour(x, y, a, [0.5])
# plt.contour(x, y, a, [0.5],colors=('r'))
plt.ylabel(r'\textbf{y}')
plt.xlabel(r'\textbf{x}')
fig1 = plt.gcf()
# plt.xlim(-0.1,5.5)
# plt.ylim(0.4,5.9)
plt.axis('off')
# fig1.set_size_inches(w=3,h=3, forward=True)

# fig1.savefig('67_sem_nothinc_exp.png', dpi=600,bbox_inches='tight', pad_inches = 0.1)
plt.show()

