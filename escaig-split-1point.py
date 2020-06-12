
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

    

def schmid_calc(b, n, T):
    global M
    n = np.dot(M, n)
    b = np.dot(M, b)
    T = T / np.linalg.norm(T)
    anglen = np.arccos(np.dot(n, T) / np.linalg.norm(n))
    angleb = np.arccos(np.dot(b, T) / np.linalg.norm(b))
    s =np.cos(anglen) * np.cos(angleb)
 
    return s


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


b=np.array([0,1,1])
n=np.array([-1,-1,1]) # Burgers vector and plane normals type in crystal coordinates

T=np.array([0,0,1]) #straining axis in fixed coordinates
phi1,phi,phi2=-110.3,64.1,76.0 #orientation
M=rotation(phi1,phi,phi2)

P=schmid_calc(b,n,T) #Schmid factor
bc=np.cross(n,b) #for positive shear (tension) the edge component of the Burgers vector is in the direction nxb, reverse for compression

diff_s=schmid_calc(bc,n,T)/P #Escaig factor t_s/sigma_a

print  bc,diff_s

