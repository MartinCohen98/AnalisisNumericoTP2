# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 14:34:29 2018

@author: Manuel
"""
#Imports
import numpy as np
import pylab as pl #Rutinas gráficas
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

xl=0; xr=1; # x domain [xl,xr]
J = 20; # J: number of division for x
dx = (xr-xl) / J; # dx: mesh size
t0 = 0.0; tf = 0.1; # initial and final simulation time
Nt = 100; # Nt: number of time steps
dt = tf/Nt;

mu = dt/(dx)**2;
#if mu > 0.5: # make sure dt satisy stability condition
#    raise ValueError('mu should < 0.5!')

# Evaluate the initial conditions
x = np.linspace(xl, xr, J+1)# generate the grid point
f = (np.sin(np.pi*x))**2

#Initialize de solution matrix in zeros
u = np.zeros([J+1,Nt+1]);

# Find the approximate solution at each time step
tt = np.linspace(t0, tf, Nt+1);
for n in np.arange(Nt+1):
    t = tt[n]; # current time
    # boundary condition at left side
    gl = 0
    # boundary condition at right side
    gr = 0
    if n==0: # first time step
        for j in np.arange(1,J): # interior nodes
            u[j,n] = f[j] + dt*f[j]*(1-f[j]) + mu*(f[j+1]-2*f[j]+f[j-1]);
        u[0,n] = gl; # the left-end point
        u[J,n] = gr; # the right-end point
    else:
        for j in np.arange(1,J): # interior nodes
            u[j,n]=u[j,n-1]+mu*(u[j+1,n-1]-2*u[j,n-1]+u[j-1,n-1]);
        u[0,n] = gl; # the left-end point
        u[J,n] = gr; # the right-end point
        
# Plot the results
print('Ploteando curvas de nivel')
pl.figure(1, figsize=(10,7))
#colormap(gray); # draw gray figure
#pl.surf(x,tt, u); # 3-D surface plot
X, Y = np.meshgrid(x, tt)
pl.contourf(X, Y, u.T, cmap=cm.coolwarm); # 3-D surface plot
pl.colorbar()
pl.xlabel('x')
pl.ylabel('t')
#pl.zlabel('u')
pl.title('Solucion numerica para u - Curvas de nivel')
pl.savefig('tp2-curvas_nivel_u.png')
pl.show()

print('Ploteando 3D')
fig = plt.figure(2, figsize=(10,7))
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, u.T, cmap=cm.coolwarm)
pl.xlabel('x')
pl.ylabel('t')
plt.title('Solucion numerica para u - Grafico 3D')
pl.savefig('tp2-grafico_3D_u.png')
plt.show()

print('Evolucion de u para distintos t')
pl.figure(3, figsize=(10,7))
#Primera y final salen en linea gruesa
#el resto en linea fina punteada
pl.plot(x, u[:, 0], 'b', label='t = ' + str(t0), lw=2)
pl.plot(x, u[:, Nt*1//4], '--', label='t = ' + str(t0 + Nt*1//4 * dt))
pl.plot(x, u[:, Nt//2], '--', label='t = ' + str(t0 + Nt//2 * dt))
pl.plot(x, u[:, Nt*3//4], '--', label='t = ' + str(t0 + Nt*3//4 * dt))
pl.plot(x, u[:, -1], 'r', label='t = ' + str(tf), lw=2)
pl.legend(loc='best')
pl.grid()
pl.xlabel('x')
pl.ylabel('u')
pl.title('Evolucion de u para distintos t')
pl.savefig('tp2-evolucion_u.png')
pl.show()