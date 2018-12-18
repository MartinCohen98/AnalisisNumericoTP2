# -*- coding: utf-8 -*-
"""
Created on Thu Dec 13 14:34:29 2018

@author: Manuel
"""
#Imports
import numpy as np
import pylab as pl #Rutinas gráficas
import matplotlib.pyplot as plt
from matplotlib import cm



def graficarEcuacionDelCalor(D, tf, J, Nt, nombre1, nombre2, nombre3) :

    xl=0; xr=1; # x domain [xl,xr]
    dx = (xr-xl) / J; # dx: mesh size
    t0 = 0.0; # initial simulation time
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
        # boundary condition at left side
        gl = 0
        # boundary condition at right side
        gr = 0
        if n==0: # first time step
            for j in np.arange(1,J): # interior nodes
                u[j,n] = f[j] + dt*f[j]*(1-f[j]) + D*mu*(f[j+1]-2*f[j]+f[j-1]);
            u[0,n] = gl; # the left-end point
            u[J,n] = gr; # the right-end point
        else:
            for j in np.arange(1,J): # interior nodes
                u[j,n]=u[j,n-1] + dt*u[j, n-1]*(1-u[j, n-1]) + D*mu*(u[j+1,n-1]-2*u[j,n-1]+u[j-1,n-1]);
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
    pl.savefig(nombre1)
    pl.show()
    
    print('Ploteando 3D')
    fig = plt.figure(2, figsize=(10,7))
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(X, Y, u.T, cmap=cm.coolwarm)
    pl.xlabel('x')
    pl.ylabel('t')
    plt.title('Solucion numerica para u - Grafico 3D')
    pl.savefig(nombre2)
    plt.show()
    
    print('Evolucion de u para distintos t')
    pl.figure(3, figsize=(10,7))
    #Primera y final salen en linea gruesa
    #el resto en linea fina punteada
    pl.plot(x, u[:, 0], 'b', label='t = ' + str(t0), lw=2)
    pl.plot(x, u[:, Nt*1//4], '--', label='t = ' + str(t0 + Nt//8 * dt))
    pl.plot(x, u[:, Nt//2], '--', label='t = ' + str(t0 + Nt//4 * dt))
    pl.plot(x, u[:, Nt*3//4], '--', label='t = ' + str(t0 + Nt*3//8 * dt))
    pl.plot(x, u[:, Nt*1//4], '--', label='t = ' + str(t0 + Nt//2 * dt))
    pl.plot(x, u[:, Nt//2], '--', label='t = ' + str(t0 + Nt*5//8 * dt))
    pl.plot(x, u[:, Nt*3//4], '--', label='t = ' + str(t0 + Nt*3//4 * dt))
    pl.plot(x, u[:, Nt//2], '--', label='t = ' + str(t0 + Nt*7//8 * dt))
    pl.plot(x, u[:, -1], 'r', label='t = ' + str(tf), lw=2)
    pl.legend(loc='best')
    pl.grid()
    pl.xlabel('x')
    pl.ylabel('u')
    pl.title('Evolucion de u para distintos t')
    pl.savefig(nombre3)
    pl.show()
    
    return pl

graficarEcuacionDelCalor(1, 0.1, 20, 100, 'CurvasDeNivelCaso1.png','SolucionNumericaCaso1.png', 'EvolucionCaso1.png')
graficarEcuacionDelCalor(1, 0.3, 15, 500, 'CurvasDeNivelCaso2.png','SolucionNumericaCaso2.png', 'EvolucionCaso2.png')
graficarEcuacionDelCalor(1, 1, 15, 500, 'CurvasDeNivelCaso3.png','SolucionNumericaCaso3.png', 'EvolucionCaso3.png')
graficarEcuacionDelCalor(1, 3, 15, 3000, 'CurvasDeNivelCaso4.png','SolucionNumericaCaso4.png', 'EvolucionCaso4.png')
graficarEcuacionDelCalor(0.5, 0.1, 20, 100, 'CurvasDeNivelCaso5.png','SolucionNumericaCaso5.png', 'EvolucionCaso5.png')
graficarEcuacionDelCalor(1.5, 0.3, 15, 500, 'CurvasDeNivelCaso6.png','SolucionNumericaCaso6.png', 'EvolucionCaso6.png')