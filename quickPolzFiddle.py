import numpy as np

pQ = 0.049594
pU = 0.027613
p = 0.056764

pQ = -0.0497
pU = -0.0045
p = 0.0499


rotateBy = 0.
# rotateBy = -54.

theta1 = 0.5 * np.arctan(pU/pQ) / np.pi * 180
theta2 = 0.5 * np.arctan2(pU, pQ) / np.pi * 180
p_calc = np.sqrt(pQ**2 + pU**2)

print('Theta 1: %f' % theta1)
print('Theta 2: %f' % theta2)
print('Calcd p: %f' % p_calc)

# newTheta = theta1 + rotateBy
newTheta = theta2 + rotateBy
pQ_new = p_calc * np.cos(2 * newTheta/180*np.pi)
pU_new = p_calc * np.sin(2 * newTheta/180*np.pi)
print('new pQ: %f' % pQ_new)
print('new pU: %f' % pU_new)
print('new theta: %f' % newTheta)