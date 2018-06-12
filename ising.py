#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
script general para hacer una corrida a un set de parámetros,
beta, tamaño de la red,
"""
import numpy as np
import matplotlib.pyplot as plt
import random as rand
plt.ion()

#############################################################

def calcMagnet(S):
	"""
	Calcula y devuelve la magnetizacion para la red de spines S
	"""
	M = np.sum(S)
	return M

def calcEnergia(S):
    H = np.zeros((L, L))
    for i in range(L):
        for j in range(L):
            k = (i + 1) % L
            l = (j + 1) % L
            H[i,j] = - S[i,j] * (S[k,j] + S[i,l])
    energia = np.sum(H)
    return energia

def difEnergia(S, i, j):
    k = (i + 1) % L
    l = (j + 1) % L
    dE = 2*S[i,j]*(S[k,j] + S[i,l] + S[i-1, j] + S[i,j-1])
    return dE

def ising2Dpaso(S, beta):
	# """
	# Calcula el proximo estado de la red de spines, para un beta dado.
	# Devluelve el próximo S, y los cambios de energia y magnetizacion
	# """
    i = rand.choice(np.arange(L))
    j = rand.choice(np.arange(L))
    dE = difEnergia(S, i, j)
    dM = -2 * S[i,j]
    if (dE <= 0 or rand.random() < np.exp(-beta*dE)):
        S[i,j] = -S[i,j]
        return S, dE, dM
    return S, 0, 0
#############################################################

#Aca defino los parámetros y corro la cadena de markov
#Lado de la red,
L = 32
T = 10000
beta = 1/T

#propongo un estado inicial al azar
#S es una matriz de 1 y -1 indicando las dos proyecciones de
#espin
S = 2*(np.random.rand(L,L)>0.5) -1;

# defino la cantidad de iteraciones de cada etapa
npre = 10000
npasos = 2000
# me genero arrays vacios, a ser llenados con los valores de energia y magnetizacion
energia= np.zeros(npasos)
magnet = np.zeros(npasos)

#pretermalizo
#ising2Dpaso hace un nuevo elemento de la cadena de Markov
#la tienen que escribir Uds...
for n in range(npre):
    S, dE, dM = ising2Dpaso(S,beta)

# muestro el estado inicial
fig_s, ax_s = plt.subplots()
ax_s.imshow(S, interpolation='none')
fig_s.show()

energia[0] = calcEnergia(S)
magnet[0] = calcMagnet(S)

# Preparamos la figura para graficar
# la energia y magnetizacion
fig_em, ax_list = plt.subplots(2,1)
ax_e, ax_m = ax_list

for n in range(npasos-1):

    S, dE, dM = ising2Dpaso(S,beta)
    energia[n+1] = energia[n] + dE
    magnet[n+1] = magnet[n] + dM

    #cada 10 pasos muestro el nuevo estado
    if n%10 == 0:

		# grafico el estado de la red, actualizandolo en cada iteracion
        ax_s.clear()
        ax_s.set_title("n=%i beta=%.2f mag=%.2f energia=%.2f"%(n,beta,magnet[n],energia[n]))
        ax_s.imshow(S, interpolation='none')
        fig_s.canvas.draw()
        plt.pause(0.01)

        # graficamos la energia y magnetizacion hasta el paso actual
        ax_e.cla()
        ax_m.cla()

        ax_e.set_xlim(0, npasos)
        ax_e.set_ylabel('Energia')
        ax_m.set_xlim(0, npasos)
        ax_m.set_ylabel('Magnetizacion')

        ax_e.grid(True)
        ax_m.grid(True)

        ax_e.plot(energia[:n])
        ax_m.plot(magnet[:n])
        fig_em.canvas.draw()
        plt.pause(0.01)
