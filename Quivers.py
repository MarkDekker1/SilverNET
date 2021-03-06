import matplotlib.pyplot as plt
import numpy as np
from numpy import ma

X, Y = np.meshgrid(np.arange(0, 2 * np.pi, .2), np.arange(0, 2 * np.pi, .2))
U = np.cos(X)
V = np.sin(Y)

U=[[1,2,3],[1,2,3],[1,2,3]]
V=[[1,2,3],[1,2,3],[1,2,3]]

A=np.meshgrid(U,V)

# 1
plt.figure()
Q = plt.quiver(A)
#qk = plt.quiverkey(Q, 0.5, 0.92, 2, r'$2 \frac{m}{s}$', labelpos='W',
#                   fontproperties={'weight': 'bold'})

l, r, b, t = plt.axis()
dx, dy = r - l, t - b
plt.axis([l - 0.05*dx, r + 0.05*dx, b - 0.05*dy, t + 0.05*dy])

plt.title('Minimal arguments, no kwargs')