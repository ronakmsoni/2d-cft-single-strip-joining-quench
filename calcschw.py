# This code takes a set of points written in output/angles.txt and calculates the Schwarzian of the doubly-connected Schwarz-Christoffel (DSC) map.
# The parameters of the map need to be manually input into the code.
#
# As written, this code only handles points at the reflection-symmetric slice in a specific DSC map (to a plane with two slits arranged in an inversion symmetric fashion), but it can be easily modified.

import math
import numpy

pi = math.pi
e = math.e

# Parameters of DSC map.
eps = .00001 #Sets the scale for the slits' ends
mu = 0.4837679155890 #annulus inner radius, from fortran. read off DSCPACK.OUTPUT
beta = [ [1,-3],  [1,1] ] #Slit exponents
pvs = [ [1,-1],[ mu, -1*mu ]  ] #Prevertices

absz = math.sqrt(mu) #pre-image of time slice
tol = 1E-20 

npt = 2000
width = 2

NConv = 25
NConv2 = 20



#Generate testing grid
nTestR = 100
nTestA = 100
testPts = []
for i in range(1,nTestR):
    rad = mu + i*(1-mu)/(nTestR+1)
    for j in range(1,nTestA):
        phi = -pi + 2*pi*j/nTestA
        testPts.append( rad * e**(phi*1j) )

#Full pre-Schwarzian f''/f'
#deLillo-Elcrat-Pfaltzgraff 2001, eqn 8
def pre_schw_full (N, z):
    preSchw = 0
    for i  in range(N): #i is nu in that equation.
        preSchw +=  1 * beta[0][0] * ( -1 * (mu**(2* i   )/pvs[0][0]) / (1 - mu**(2* i   ) * z/pvs[0][0]) + ( mu**(2*(i+1)) * pvs[0][0]/z**2 ) / ( 1 - mu**(2*(i+1)) * pvs[0][0]/z ) ) 
        preSchw +=  1 * beta[0][1] * ( -1 * (mu**(2* i   )/pvs[0][1]) / (1 - mu**(2* i   ) * z/pvs[0][1]) + ( mu**(2*(i+1)) * pvs[0][1]/z**2 ) / ( 1 - mu**(2*(i+1)) * pvs[0][1]/z ) ) 
        preSchw +=      beta[1][0] * ( -1 * (mu**(2*(i+1))/pvs[1][0]) / (1 - mu**(2*(i+1)) * z/pvs[1][0]) + ( mu**(2* i   ) * pvs[1][0]/z**2 ) / ( 1 - mu**(2* i   ) * pvs[1][0]/z ) ) 
        preSchw +=      beta[1][1] * ( -1 * (mu**(2*(i+1))/pvs[1][1]) / (1 - mu**(2*(i+1)) * z/pvs[1][1]) + ( mu**(2* i   ) * pvs[1][1]/z**2 ) / ( 1 - mu**(2* i   ) * pvs[1][1]/z ) ) 
    return preSchw
    
#error check
#errors = []
#for pt in testPts:
#    if abs(pre_schw_full(NConv,pt) - pre_schw_full(NConv+1,pt)) > tol:
#        errors.append(abs(pre_schw_full(NConv,pt) - pre_schw_full(NConv+1,pt)))
#print(errors)

#Calculates the Schwarzian
def schw(N,z):
    sch = - .5 * pre_schw_full(N,z)**2
    for i  in range(N):
        sch +=  1 * beta[0][0] * ( -1 * (mu**(4* i   )/pvs[0][0]**2) / (1 - mu**(2* i   ) * z/pvs[0][0])**2 - mu**(2*(i+1)) * pvs[0][0]/z**3 * (2 - mu**(2*(i+1)) * pvs[0][0]/z ) / ( 1 - mu**(2*(i+1)) * pvs[0][0]/z )**2 ) 
        sch +=  1 * beta[0][1] * ( -1 * (mu**(4* i   )/pvs[0][1]**2) / (1 - mu**(2* i   ) * z/pvs[0][1])**2 - mu**(2*(i+1)) * pvs[0][1]/z**3 * (2 - mu**(2*(i+1)) * pvs[0][1]/z ) / ( 1 - mu**(2*(i+1)) * pvs[0][1]/z )**2 ) 
        sch +=  1 * beta[1][0] * ( -1 * (mu**(4*(i+1))/pvs[1][0]**2) / (1 - mu**(2*(i+1)) * z/pvs[1][0])**2 - mu**(2* i   ) * pvs[1][0]/z**3 * (2 - mu**(2* i   ) * pvs[1][0]/z ) / ( 1 - mu**(2* i   ) * pvs[1][0]/z )**2 ) 
        sch +=  1 * beta[1][1] * ( -1 * (mu**(4*(i+1))/pvs[1][1]**2) / (1 - mu**(2*(i+1)) * z/pvs[1][1])**2 - mu**(2* i   ) * pvs[1][1]/z**3 * (2 - mu**(2* i   ) * pvs[1][1]/z ) / ( 1 - mu**(2* i   ) * pvs[1][1]/z )**2 ) 
    return sch

#error check
#errors = []
#for pt in testPts:
#    if abs(schw(NConv,pt) - schw(NConv+1,pt)) > tol:
#        errors.append(abs(schw(NConv,pt) - schw(NConv+1,pt)))
#print(errors)
#babadook

#for pt in testPts:
#    print schw(NConv,pt)

f = open("output/angles.txt","r")
angs = f.read().splitlines()
f.close()

for i in range(len(angs)):
    angs[i] = float(angs[i])

#Actual calculation
f = open("output/schws.txt","w")
f2 = open("output/pre-schws.txt","w")
f3 = open("output/energy-terms.txt","w")
for ang in angs:
    z = absz * e**(1j*ang)
    sch = schw(NConv,z)
    f.write(str(sch.real))
    f.write(" + I * ")
    f.write(str(sch.imag))
    f.write("\n")
    f2.write(str(pre_schw_full(NConv,z)))
    f2.write("\n")
    f3.write(str((z*z*sch).real))
    f3.write(" + I * ")
    f3.write(str((z*z*sch).imag))
    f3.write("\n")

f.close()
f2.close()
f3.close()
