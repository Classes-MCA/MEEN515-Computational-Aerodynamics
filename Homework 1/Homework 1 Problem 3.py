# This program will solve the Leapfrogging Vortex Problem

# Each vortex induces a velocity on each other vortex as Vind = Gamma/(2*pi*distance)
# In vector terms: Vind = Gamma x r / (2*pi*||r||^2)

import numpy as np
import scipy.linalg as la
import matplotlib.pyplot as plt

animate = False;

vortex_A = np.array([0,0,0]) # Bottom Left
vortex_B = np.array([0,1,0]) # Top Left
vortex_C = np.array([1,1,0]) # Top Right
vortex_D = np.array([1,0,0]) # Bottom Right

# For efficiency, we will store the Gamma value in the third dimension of the 
# vortices. Remember that indices are zero-based, so '2' means the third index
vortex_A[2] = -1
vortex_B[2] = 1
vortex_C[2] = 1
vortex_D[2] = -1

startTime = 0
endTime = 50
timeSteps = (endTime - startTime) * 100

t,dt = np.linspace(startTime,endTime,timeSteps,retstep = True)

def inducedVelocity(currentVortex,otherVortex1,otherVortex2,otherVortex3):

    otherVortex1_Strength = [0,0,otherVortex1[2]]
    otherVortex2_Strength = [0,0,otherVortex2[2]]
    otherVortex3_Strength = [0,0,otherVortex3[2]]
    
    inducedVelocity1 = np.cross(otherVortex1_Strength,currentVortex - otherVortex1)/(la.norm(currentVortex - otherVortex1)**2)
    inducedVelocity2 = np.cross(otherVortex2_Strength,currentVortex - otherVortex2)/(la.norm(currentVortex - otherVortex2)**2)
    inducedVelocity3 = np.cross(otherVortex3_Strength,currentVortex - otherVortex3)/(la.norm(currentVortex - otherVortex3)**2)
    
    inducedVelocity = inducedVelocity1 + inducedVelocity2 + inducedVelocity3
    inducedVelocity = inducedVelocity/(2*np.pi)
    
    return inducedVelocity

vortex_A_History = np.zeros((len(t),3))
vortex_B_History = np.zeros((len(t),3))
vortex_C_History = np.zeros((len(t),3))
vortex_D_History = np.zeros((len(t),3))

for i in range(len(t)): 
    
    # Calculate the induced velocities at each vortex due to the other vortices
    V_A = inducedVelocity(vortex_A,vortex_B,vortex_C,vortex_D)
    V_B = inducedVelocity(vortex_B,vortex_A,vortex_C,vortex_D)
    V_C = inducedVelocity(vortex_C,vortex_B,vortex_A,vortex_D)
    V_D = inducedVelocity(vortex_D,vortex_B,vortex_C,vortex_A)
    
    # Calculate the displacement of each vortex during this time step
    vortex_A = vortex_A + V_A*dt
    vortex_B = vortex_B + V_B*dt
    vortex_C = vortex_C + V_C*dt
    vortex_D = vortex_D + V_D*dt
    
    #currentPoints = np.array([vortex_A[0], vortexA_[1]; vortex_B[0], vortex_B[1]])
    
    if animate == True:
        plt.clf()
        plt.scatter(vortex_A[0],vortex_A[1],color = 'blue')
        plt.scatter(vortex_B[0],vortex_B[1],color = 'blue')
        plt.scatter(vortex_C[0],vortex_C[1],color = 'red')
        plt.scatter(vortex_D[0],vortex_D[1],color = 'red')
        plt.ylim([-0.5,1.5])
        plt.xlim([-0.5,3])
        plt.title('Leapfrogging Vortex Rings')
        plt.xlabel('x (length units)')
        plt.ylabel('y (length units)')
        plt.draw()
        plt.pause(0.001)
        
    vortex_A_History[i] = vortex_A
    vortex_B_History[i] = vortex_B
    vortex_C_History[i] = vortex_C
    vortex_D_History[i] = vortex_D
    
plt.figure()
plt.plot(vortex_A_History[:,0],vortex_A_History[:,1],color = 'blue')
plt.plot(vortex_B_History[:,0],vortex_B_History[:,1],color = 'blue')
plt.plot(vortex_C_History[:,0],vortex_C_History[:,1],color = 'red')
plt.plot(vortex_D_History[:,0],vortex_D_History[:,1],color = 'red')
plt.title('Leapfrogging Vortex Rings')
plt.xlabel('x (length units)')
plt.ylabel('y (length units)')