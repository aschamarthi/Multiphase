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
plt.rcParams.update({'font.size': 12})
plt.rcParams['axes.xmargin'] = 0


# 247 \mu seconds

x,y,d,u,v,p,a,grad,phi=np.loadtxt('Rslt0026.plt', delimiter=None, unpack=True,skiprows=3)

# x,y,d,u,v,p,a,grad,phi=np.loadtxt('Rslt0035.plt', delimiter=None, unpack=True,skiprows=3)

g=512*4
k=768*4
x = x.reshape(g,k)
y = y.reshape(g,k)
# d = d.reshape(g,k)
# p = p.reshape(g,k)
# u = u.reshape(g,k)
# v = v.reshape(g,k)
# a = a.reshape(g,k)
# phi = phi.reshape(g,k)
grad = grad.reshape(g,k)


# gop=np.exp(-20000*grad/grad.max())

# blah=np.log(abs(grad)+1)

# t=p/(d*287)

# print(t.max())

# plt.contourf(x,y,p,40,cmap='jet')



plt.imshow(grad, vmin = 1, vmax = 2400, cmap=plt.cm.Blues, origin='upper', 
           extent=[x.min(), x.max(), y.min(), y.max()])
# plt.contour(x,y,a,[0.99],colors='r')
# plt.imshow(blah, vmin = 0, vmax = 15, cmap=plt.cm.gray_r, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()])
# plt.contour(x,y,a,1,cmap='gray',linewidth=0.5)
# plt.imshow(gop, vmin = 0, vmax = 1.0, cmap=plt.cm.Blues_r, origin='lower', 
           # extent=[x.min(), x.max(), y.min(), y.max()])
# plt.contour(x, y, a, [0.5],colors=('r'),linewidths=0.5)
plt.arrow(0.051, 0.049, 0.01725, 0.00, fc="k", ec="r",head_width=0.002, head_length=0.002, width=.0002)

plt.ylabel(r'\textbf{y}')
plt.xlabel(r'\textbf{x}')
fig1 = plt.gcf()

plt.axis('off')
# plt.xlim(0.0,0.08)
# fig1.set_size_inches(w=8,h=5, forward=True)
# plt.colorbar(orientation="vertical")
fig1.savefig('chx67_arrow.png', dpi=600,bbox_inches='tight', pad_inches = 0.0)
plt.show()

# 318 \mu seconds

