
from __future__ import division
import numpy as np
from matplotlib import pyplot as plt


def proj(x, y, z):

    if z == 1:
        X = 0
        Y = 0
    elif z < -0.000001:
        X = -x / (1 - z)
        Y = -y / (1 - z)
    else:

        X = x / (1 + z)
        Y = y / (1 + z)

    return np.array([X, Y], float)


############################################
# Rotation Euler
#
###########################################

def rotation(phi1, phi, phi2):
    phi1 = phi1 * np.pi / 180
    phi = phi * np.pi / 180
    phi2 = phi2 * np.pi / 180
    R = np.array([[np.cos(phi1) * np.cos(phi2) - np.cos(phi) * np.sin(phi1) * np.sin(phi2),
                   -np.cos(phi) * np.cos(phi2) * np.sin(phi1) - np.cos(phi1) *
                   np.sin(phi2), np.sin(phi) * np.sin(phi1)], [np.cos(phi2) * np.sin(phi1)
                                                               + np.cos(phi) * np.cos(phi1) * np.sin(phi2), np.cos(phi) * np.cos(phi1)
                                                               * np.cos(phi2) - np.sin(phi1) * np.sin(phi2), -np.cos(phi1) * np.sin(phi)],
                  [np.sin(phi) * np.sin(phi2), np.cos(phi2) * np.sin(phi), np.cos(phi)]], float)
    return R

###################################################################
#
# Rotation around a given axis
#
##################################################################


def Rot(th, a, b, c):
    th = th * np.pi / 180
    no = np.linalg.norm([a, b, c])
    aa = a / no
    bb = b / no
    cc = c / no
    c1 = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    c2 = np.array([[aa**2, aa * bb, aa * cc], [bb * aa, bb**2, bb * cc], [cc * aa,
                                                                          cc * bb, cc**2]], float)
    c3 = np.array([[0, -cc, bb], [cc, 0, -aa], [-bb, aa, 0]], float)
    R = np.cos(th) * c1 + (1 - np.cos(th)) * c2 + np.sin(th) * c3

    return R
    

def schmid_calc(b, n, T):
    global M
    n = np.dot(M, n)
    b = np.dot(M, b)
    T = T / np.linalg.norm(T)
    anglen = np.arccos(np.dot(n, T) / np.linalg.norm(n))
    angleb = np.arccos(np.dot(b, T) / np.linalg.norm(b))
    s =np.cos(anglen) * np.cos(angleb)

    return s

####################
#
# Get Shockley partials. Use the fact that for a given normal the leading partial is the one closest to the <001> direction between n and b. Equivalent to say that the leading partial is the one toward -n x xi in the Tompson tetrahedron.
#
#########################
def shockley(b,n):
	if b[0]==0:
		c1=np.array([0,0,1])
		c2=np.array([0,1,0])
	
	if b[1]==0:
		c1=np.array([1,0,0])
		c2=np.array([0,0,1])
	
	if b[2]==0:
		c1=np.array([1,0,0])
		c2=np.array([0,1,0])	
	
	if np.dot(c1,b)>0 and np.dot(c1,n)>0:
			r=np.cross(n,c1)
	elif np.dot(c1,b)<0 and np.dot(c1,n)<0:
			r=np.cross(n,-c1)
	else:
		if np.dot(c2,b)>0 and np.dot(c2,n)>0:
			r=np.cross(n,c2)
		elif np.dot(c2,b)<0 and np.dot(c2,n)<0:
			r=np.cross(n,-c2)
	
	sl=np.dot(Rot(90,r[0],r[1],r[2]),n)*np.sqrt(2)
	st=6*(b/2-sl/6)	
	
	return st,sl


def unique_rows(a):
    a = np.ascontiguousarray(a)
    unique_a = np.unique(a.view([('', a.dtype)] * a.shape[1]))
    return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

#################"
#
# Get equivalent pole
#
#################

def schmid_pole(n):
	pole1,pole2,pole3=n[0],n[1],n[2]
	N = np.array([pole1, pole2, pole3])
	N = np.vstack((N, np.array([pole1, pole2, -pole3])))
	N = np.vstack((N, np.array([pole1, -pole2, pole3])))
	N = np.vstack((N, np.array([-pole1, pole2, pole3])))
	N = np.vstack((N, np.array([pole2, pole1, pole3])))
	N = np.vstack((N, np.array([pole2, pole1, -pole3])))
	N = np.vstack((N, np.array([pole2, -pole1, pole3])))
	N = np.vstack((N, np.array([-pole2, pole1, pole3])))
	N = np.vstack((N, np.array([pole2, pole3, pole1])))
	N = np.vstack((N, np.array([pole2, pole3, -pole1])))
	N = np.vstack((N, np.array([pole2, -pole3, pole1])))
	N = np.vstack((N, np.array([-pole2, pole3, pole1])))
	N = np.vstack((N, np.array([pole1, pole3, pole2])))
	N = np.vstack((N, np.array([pole1, pole3, -pole2])))
	N = np.vstack((N, np.array([pole1, -pole3, pole2])))
	N = np.vstack((N, np.array([-pole1, pole3, pole2])))
	N = np.vstack((N, np.array([pole3, pole1, pole2])))
	N = np.vstack((N, np.array([pole3, pole1, -pole2])))
	N = np.vstack((N, np.array([pole3, -pole1, pole2])))
	N = np.vstack((N, np.array([-pole3, pole1, pole2])))
	N = np.vstack((N, np.array([pole3, pole2, pole1])))
	N = np.vstack((N, np.array([pole3, pole2, -pole1])))
	N = np.vstack((N, np.array([pole3, -pole2, pole1])))
	N = np.vstack((N, np.array([-pole3, pole2, pole1])))
	return N

###############
#
# Get all the SF and remove duplicate
#
####################
def schmid(B, N, T):
    
    P = np.array([0, 0, 0, 0, 0, 0, 0])
    for i in range(0, np.shape(N)[0]):
        for j in range(0, np.shape(B)[0]):
        	if np.abs(np.dot(N[i, :], B[j, :])) < 0.0001:
		        s = schmid_calc(B[j, :], N[i, :], T)
		        R = np.array([s, N[i, 0], N[i, 1], N[i, 2], B[j, 0], B[j, 1], B[j, 2]])
		        P = np.vstack((P, R))

    P = np.delete(P, (0), axis=0)
    P = unique_rows(P)
    return P

############################
#
# Main
#
##############################"

f= open("escaig-result-compression.txt","w+") #file where results will be written

b=np.array([1,1,0])
n=np.array([1,1,1]) #primary Burgers and plane normals type in crystal coordinates
b=schmid_pole(b)
n=schmid_pole(n)
T=np.array([0,0,1]) #straining axis in fixed coordinates
phi1=0

#####################"
#
# Ranging in phi1 and phi2, which describes a little larger standard triangle, nn to cover a larger area , nm to increase resolution
#
#####################"
nn=1
nm=1
for phi in range(0,56*nn):
	for phi2 in range(0,46*nn):
		M=rotation(phi1,phi/nm,phi2/nm)
		if np.dot(M.T[:,2],[0,-1,1])>-0.006: #to get the standard triangle
			P=schmid(b,n,T)
			Pmax=P[np.argmax(np.abs(P[:,0])),:] #get primary slip system with highest SF
			bmax=Pmax[4:]
			nmax=Pmax[1:4]
			CS=np.delete(P,np.argmax(np.abs(P[:,0])),axis=0) #get the cross slip plane
			CS=CS[(CS[:,4:]==bmax).all(axis=1)][0]
			ncs=CS[1:4]
			sl,st=shockley(bmax,nmax)
			s1=schmid_calc(sl,nmax,T)
			s2=schmid_calc(st,nmax,T)
			diff_s=s1-s2
			slcs,stcs=shockley(bmax,ncs)
			s1cs=schmid_calc(slcs,ncs,T)
			s2cs=schmid_calc(stcs,ncs,T)
			diff_scs=s1cs-s2cs
			Z=M.T[:,2] #get the straining axis coordinate in crystal frame
			Zp=proj(Z[0],Z[1],Z[2]) #stereographic projection
			
			f.write(str(np.abs(Pmax[0]))+','+str(diff_s)+','+str(diff_scs)+','+str(Zp[0])+','+str(Zp[1])+','+str(phi1)+','+str(phi/nn)+','+str(phi2/nn)+'\n')
	

f.close()

