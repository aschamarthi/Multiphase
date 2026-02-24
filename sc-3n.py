# from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import colors
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
plt.rcParams.update({'font.size': 12})



# x,y,phi,vorticity=np.loadtxt('new_triple_migi_50.plt', delimiter=None, unpack=True,skiprows=2)
x,y,sig_x=np.loadtxt('sen3.plt', delimiter=None, unpack=True,skiprows=2)
# x,y,d,u,v,p,a,grad,phi=np.loadtxt('Rslt0026.plt', delimiter=None, unpack=True,skiprows=3)

g=512*4
k=768*4
# 

x = x.reshape(g,k)
y = y.reshape(g,k)
# d = d.reshape(g,k)
# p = p.reshape(g,k)
# u = u.reshape(g,k)
# v = v.reshape(g,k)
# a = a.reshape(g,k)
# vorticity = vorticity.reshape(g,k)
sig_x = sig_x.reshape(g,k)
# sig_y = sig_y.reshape(g,k)
# grad = grad.reshape(g,k)
# phi = phi.reshape(g,k)
# blah=np.log(abs(grad)+1)
# print(a.max())

# print(vorticity.min())

# plt.contourf(x,y,phi,40,cmap='jet_r')
# plt.contourf(x,y,a,40,cmap='jet_r')
plt.contour(x,y,sig_x,15,linewidths=0.25,colors=('b'))
# plt.contour(x,y,sig_y,40,linewidths=0.25,colors=('r'))
# plt.contour(x,y,u,32,linewidths=0.5,colors=('k'))
# plt.imshow(vorticity, vmin = 0.8, vmax = 7, cmap=plt.cm.Blues, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()])
# plt.imshow(blah, vmin = 0, vmax = 6, cmap=plt.cm.gray_r, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()])

# plt.contourf(vorticity, cmap=plt.cm.Spectral_r, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()],norm = colors.LogNorm(vmin = 0.8, vmax = 10))
# plt.imshow(a, vmin = a.min(), vmax = a.max(), cmap=plt.cm.Blues, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()],alpha=0.95)

# t=1-a

# plt.imshow(p, vmin = p.min(), vmax = p.max(), cmap=plt.cm.Greens, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()],alpha=0.45)

# flop = np.exp(-10*abs(grad)/grad.max())

# print(grad.max())
# print(a.max())
# print(a.min())
# plt.imshow(grad, vmin = 1, vmax = 3400, cmap=plt.cm.Blues, origin='upper', 
           # extent=[x.min(), x.max(), y.min(), y.max()])
# print(blah.max())
# plt.imshow(flop, vmin = 0.0, vmax = 1.0, cmap=plt.cm.Blues_r, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()])
# plt.imshow(a, vmin = 0, vmax = 1, cmap=plt.cm.jet_r, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()])

# plt.imshow(phi, vmin = phi.min(), vmax = phi.max(), cmap=plt.cm.gray_r, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()])
plt.ylabel(r'\textbf{y}')
plt.xlabel(r'\textbf{x}')
fig1 = plt.gcf()
# plt.xlim(2.3,7.0)
# plt.ylim(0.0,3.0)

# plt.arrow(3.25, 0.5, -0.5725, 0.0, fc="k", ec="r",head_width=0.03, head_length=0.05)
# plt.minorticks_on()
# plt.grid(b=True, which='minor')
# plt.grid(b=True, which='major')



plt.axis('off')
# plt.spines['right'].set_visible(False)
# plt.spines['top'].set_visible(False)
# fig1.set_size_inches(w=4,h=4, forward=True)

fig1.savefig('noise.png', dpi=300,bbox_inches='tight', pad_inches = 0.0)
# fig1.savefig('trip_mig_50_vorticity.pdf', dpi=1200,bbox_inches='tight', pad_inches = 0.1)
plt.show()

