
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
# Get Shockley partials. 
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
	
	st=np.dot(Rot(90,r[0],r[1],r[2]),n)*np.sqrt(2)
	sl=6*(b/2-st/6)	
	
	return sl,st

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


b=np.array([0,-1,1])
n=np.array([1,-1,-1]) #primary Burgers and plane normals type in crystal coordinates

T=np.array([0,1,0]) #straining axis in fixed coordinates
phi1,phi,phi2=-179.0,62.3,118.7

M=rotation(phi1,phi,phi2)


SF=schmid_calc(b,n,T)
sl,st=shockley(b,n)
s1=schmid_calc(sl,n,T)
s2=schmid_calc(st,n,T)
diff_s=s1-s2

print 'Schmid factor: ',SF
print 'Escaig factor:', diff_s
print 'partials (leading and trailing): ',sl,st
print 'Schmid factor on leading and trailing: ', s1,s2


