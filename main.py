# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 09:53:20 2020

@author: marco
"""
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
#import configparser
#import sys
#from sys import argv
import functions


config = configparser.ConfigParser()
config.read('input.txt')
#config.read(sys.argv[1])

N = int( config.get('parameters', 'N') )
iterations = int( config.get('parameters', 'iterations') )
nstep = int( config.get('parameters', 'nstep') )
T_min = float( config.get('parameters', 'T_min') )
T = float( config.get('parameters', 'T') )
alpha = float( config.get('parameters', 'alpha') )

'''
city = np.zeros((2,N))
citynum = list(range(N))
Tem = []
accelist = []
for i in range(N):
    for j in (0,1):
        city[j][i] = rand.rand() #cities coordinates 
'''

'''
OPTIONS (RAW)
Moves: swap, block reverse, prune and graft
Criteria: Metropolis, distance
'''
iniziale=functions.lenght(N)
print('Distanza iniziale:',iniziale)
corrente=[iniziale]
functions.plpath(N)

for ii in range(iter):
    
    functions.anneal_BRev_distance(N, alpha, T)
    #functions.anneal_BRev_Metropolis(N, alpha, T)
    #functions.anneal_swap_Metropolis()
    #functions.anneal_swap_distance(N, alpha, T)
    #functions.anneal_PG_Metropolis(N, alpha, T)
    #functions.anneal_PG_distance(N, alpha, T)
    
    if ii==0:
        functions.confronto()
    corrente.append(functions.lenght(N))
    print('Distanza corrente:',corrente[ii+1],', Iter:',ii+1)
    functions.plpath(N)
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