#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import random as rand
plt.ion()

def calcMagnet(S):
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
	i = rand.choice(np.arange(L))
	j = rand.choice(np.arange(L))
	r = rand.random()
	dE = difEnergia(S, i, j)
	dM = -2 * S[i,j]
	if dE <= 0:
		S[i,j] = -S[i,j]
		return S, dE, dM
	elif r < np.exp(-beta*dE):
		S[i,j] = -S[i,j]
		return S, dE, dM
	return S, 0, 0

L = 32
T = 3
beta = 1.0/T

S = 2*(np.random.rand(L,L)>0.5) -1;

npre = 1000000
npasos = 1000

energia= np.zeros(npasos)
magnet = np.zeros(npasos)

for n in range(npre):
    S, dE, dM = ising2Dpaso(S,beta)

# fig_s, ax_s = plt.subplots()
# ax_s.imshow(S, interpolation='none')
# fig_s.show()

energia[0] = calcEnergia(S)
magnet[0] = calcMagnet(S)

# fig_em, ax_list = plt.subplots(2,1)
# ax_e, ax_m = ax_list

for n in range(npasos-1):

    S, dE, dM = ising2Dpaso(S,beta)
    energia[n+1] = energia[n] + dE
    magnet[n+1] = magnet[n] + dM

    # if n%10 == 0:
	#
    #     ax_s.clear()
    #     ax_s.set_title("n=%i beta=%.2f mag=%.2f energia=%.2f"%(n,beta,magnet[n],energia[n]))
    #     ax_s.imshow(S, interpolation='none')
    #     fig_s.canvas.draw()
    #     plt.pause(0.01)
	#
    #     ax_e.cla()
    #     ax_m.cla()
	#
    #     ax_e.set_xlim(0, npasos)
    #     ax_e.set_ylabel('Energia')
    #     ax_m.set_xlim(0, npasos)
    #     ax_m.set_ylabel('Magnetizacion')
	#
    #     ax_e.grid(True)
    #     ax_m.grid(True)
	#
    #     ax_e.plot(energia[:n])
    #     ax_m.plot(magnet[:n])
    #     fig_em.canvas.draw()
    #     plt.pause(0.01)
