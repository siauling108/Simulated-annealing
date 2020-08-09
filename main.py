# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 09:53:20 2020

@author: marco
"""
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
import configparser
import sys
from sys import argv
import functions


config = configparser.ConfigParser()
config.read(sys.argv[1])

N = config.get('parameters', 'N')
iterations = config.get('parameters', 'iterations')
T_min = config.get('parameters', 'T_min')
T = config.get('parameters', 'T')

city = np.zeros((2,N))
lista = list(range(N))
#lista=list(lista)
Tem = []
accelist = []
for i in range(N):
    for j in (0,1):
        city[j][i] = rand.rand() #cities coordinates 

'''
OPTIONS (RAW)
Moves: swap, block reverse, prune and graft
Criteria: Metropolis, distance
'''
iniziale=functions.lung()
print('Distanza iniziale:',iniziale)
corrente=[iniziale]
functions.plpath()
pat=[0] #andranno appese tute le lunghezze per ogni iterazione, non mi pare l abbia utilizzato

for ii in range(iter):
    
    functions.anneal_BRev_distance()
    #functions.anneal_BRev_Metropolis()
    #functions.anneal_swap_Metropolis()
    #functions.anneal_swap_distance()
    #functions.anneal_PG_Metropolis()
    #functions.anneal_PG_distance()
    
    if ii==0:
        functions.confronto()
    corrente.append(functions.lung())
    print('Distanza corrente:',corrente[ii+1],', Iter:',ii+1)
    functions.plpath()
if iter>1:
    plt.figure(2)
    plt.title("Distanza corrente in funzione del numero di iterazioni")
    plt.xlabel("Iter")
    plt.ylabel("Distanza corrente (unita' arbitrarie)")
    plt.plot(corrente,label="Distanza")
    plt.legend()
    plt.grid()
print('Distanza finale:',corrente[-1])
plt.show()
print("Diminuzione percentuale della distanza", (min(corrente)/max(corrente))*100., "%")