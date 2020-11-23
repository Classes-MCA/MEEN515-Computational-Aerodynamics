import matplotlib.pyplot as plt
import numpy as np

baseSize = np.array([100,90,50,20,10,4,2,1,0.5,0.2,0.1,0.05])

liftCoefficient = np.array([84.666,85.796,92.215,95.066,96.821,99.252,99.57,100.035,100.27,100.548,100.612,100.858])

dragCoefficient = np.array([5.479,4.679,2.891,1.449,1.085,0.569,0.383,0.337,0.285,0.264,0.251,0.249])

plt.figure()
plt.plot(1/baseSize,liftCoefficient,'b*')
plt.title('Lift Coefficient')
plt.xlabel('Base Size$^{-1}$')
plt.ylabel('Lift Coefficient')
plt.xscale('log')
plt.grid('on')

plt.figure()
plt.plot(1/baseSize,dragCoefficient,'b*')
plt.title('Drag Coefficient')
plt.xlabel('Base Size$^{-1}$')
plt.ylabel('Drag Coefficient')
plt.xscale('log')
plt.grid('on')