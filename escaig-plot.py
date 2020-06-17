from __future__ import division
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors
import matplotlib.cm as cm

    
data=np.genfromtxt('escaig-result.txt',delimiter=',') #get data from file
figure=plt.figure()
a=figure.add_subplot('111')

###########################
#Comment/uncomment to plot (for td and td', the 0 corresponds to the white lines)
############################
#######################
# t primary

#sc=a.tricontourf(data[:,3],data[:,4], data[:,0], 101,cmap=plt.cm.seismic)

###################
# td cross slip 

#r=np.max((np.abs(data[:,2].max()),np.abs(data[:,2].min())))
#levels=np.linspace(-r, r,num=255)
#sc=a.tricontourf(data[:,3],data[:,4], -data[:,2], levels=levels,cmap=plt.cm.seismic)

###################
# td' primary 

#r=np.max((np.abs(data[:,1].max()),np.abs(data[:,1].min())))
#levels=np.linspace(-r, r,num=255)
#sc=a.tricontourf(data[:,3],data[:,4], data[:,1], levels=levels,cmap=plt.cm.seismic)


######################
#orientation factor
m=((-2*data[:,2]/3+data[:,1])/data[:,0])
sc=a.tricontourf(data[:,3],data[:,4], m,301,cmap=plt.cm.seismic)
#a.tricontour(data[:,3],data[:,4], m, colors='k')




plt.colorbar(sc)
a.axis('equal')
a.axis('off')
plt.show()
