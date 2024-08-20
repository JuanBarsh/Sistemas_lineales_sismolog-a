# -*- coding: utf-8 -*-
"""
Created on Tue Apr 16 07:41:16 2024

@author: juanb
"""

# # Limpiar todas las variables
# import gc
# gc.collect()

import matplotlib.pyplot as plt
# plt.close('all') #Cerrar todas la figuras

# # Limpiar la consola
# import os
# os.system('cls')


import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy import stats
from scipy import optimize as opt
from scipy import fftpack

dt = 0.01;

starl = 110 #Desde esta línea empieza a leer
with open('BA4920170919181440','r') as acele:
    acelerograma = acele.readlines()[starl-1:]
    for line in acelerograma:
        print(line)
        
# Dividir cada string de la variable tipo list
acel_div = [cadena.split() for cadena in acelerograma]

matriz_acel = np.array([[float(valor.strip()) for valor in fila] for fila in acel_div]);

print('Tenemos la matriz construida')

AcelNS = []
AcelV = []
AcelEW = []

for i in range(0,len(matriz_acel)):
    AcelNS.append(matriz_acel[i][0])
    AcelV.append(matriz_acel[i][1])
    AcelEW.append(matriz_acel[i][2])

#Graficar los componentes
tempo = []
t = 0
dt = 0.01
for i in range(0,len(matriz_acel)):
    tempo.append([t + dt * i])

plt.show(1)
plt.plot(tempo,AcelNS)
plt.xlabel('Tiempo(s)')
plt.ylabel('Aceleración(cm/s^2)')
plt.title('Aceleración componente N-S')
plt.savefig ('NS.png', dpi=200, format='png')

plt.show(2)
plt.plot(tempo,AcelV)
plt.xlabel('Tiempo(s)')
plt.ylabel('Aceleración(cm/s^2)')
plt.title('Aceleración componente Vertical')
plt.savefig ('Vert.png', dpi=200, format='png')

plt.show(3)
plt.plot(tempo,AcelEW)
plt.xlabel('Tiempo(s)')
plt.ylabel('Aceleración(cm/s^2)')
plt.title('Aceleración componente E-W')
plt.savefig ('EW.png', dpi=200, format='png')

N = len(AcelNS);
FF = 1/(N*dt);
FNy = 1/(2*dt); #Frecuencia de Nyquist
F = np.arange(-FNy,FNy,FF);

####    Procedemos a calcular la transformada de Fourier para cada una de nuestras componentes
FAcelNS = fftpack.fft(AcelNS);

FAcelV = fftpack.fft(AcelV);

FAcelEW = fftpack.fft(AcelEW);

G0 = np.fft.fftshift(FAcelNS);
G1 = np.fft.fftshift(FAcelV);
G2 = np.fft.fftshift(FAcelEW);
plt.show(4)
plt.loglog(F,abs(G0*2/N))
plt.xlabel('Frecuencia')
plt.ylabel('Amplitud')
plt.title('Espectro de Fourier NS')
plt.savefig ('EspectroNS.png', dpi=200, format='png')

plt.show(5)
plt.loglog(F,abs(G1*2/N))
plt.xlabel('Frecuencia')
plt.ylabel('Amplitud')
plt.title('Espectro de Fourier NS')
plt.savefig ('EspectroV.png', dpi=200, format='png')

plt.show(6)
plt.loglog(F,abs(G2*2/N))
plt.xlabel('Frecuencia')
plt.ylabel('Amplitud')
plt.title('Espectro de Fourier NS')
plt.savefig ('EspectroEW.png', dpi=200, format='png')