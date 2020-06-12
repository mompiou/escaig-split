from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors
import matplotlib.cm as cm

    
data=np.genfromtxt('escaig-result.txt',delimiter=',') #get data from file
figure=plt.figure()
a=figure.add_subplot('111')

m=data[:,1]
r=np.max((np.abs(m.max()),np.abs(m.min())))
levels=np.linspace(-r, r,num=255)
sc=a.tricontourf(data[:,2],data[:,3], m, levels=levels,cmap=plt.cm.seismic)


plt.colorbar(sc)
a.axis('equal')
a.axis('off')
plt.show()
