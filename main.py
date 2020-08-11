# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 09:53:20 2020

@author: marco
"""
import numpy as np
import numpy.random as rand
import matplotlib.pyplot as plt
import configparser
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

Tem = []
accelist = []
city = functions.travel(N)
citynum = list(range(N))
start=functions.lenght(N, city)
print('Initial distance:',start)
distances=[start]
functions.plpath(N, city)

'''
OPTIONS (RAW)
Moves: swap, block reverse, prune and graft
Criteria: Metropolis, distance
'''

for ii in range(iterations):
    functions.anneal_BRev_distance(N, alpha, T, city, T_min, nstep, Tem, accelist)
    #functions.anneal_BRev_Metropolis(N, alpha, T, city, T_min, nstep, Tem, accelist) 
    #functions.anneal_swap_Metropolis(N, alpha, T, city, T_min, nstep, Tem, accelist) #comp. demanding
    #functions.anneal_swap_distance(N, alpha, T, city, T_min, nstep, Tem, accelist) #ratio vs T??
    #functions.anneal_PG_Metropolis(N, alpha, T, city, T_min, nstep, Tem, accelist, citynum) # NON FUNZ
    #functions.anneal_PG_distance(N, alpha, T, city, T_min, nstep, Tem, accelist, citynum) # NON FUNZ
    
    if ii==0:
        functions.confronto(accelist, Tem)
    distances.append(functions.lenght(N, city))
    print('Current distance:',distances[ii+1],', Iteration:',ii+1)
    functions.plpath(N, city)
if iterations>1:
    plt.figure(2)
    plt.title("Current distance vs. # of iterations")
    plt.xlabel("Iteration")
    plt.ylabel("Current distance (arb. units)")
    plt.plot(distances)
    plt.legend()
    plt.grid()
print('Final distance:',distances[-1])
plt.show()
print("Percentage decrease of the distance", (max(distances)-min(distances))/max(distances)*100., "%")
#print("Percentage decrease of the distance", (min(distances)/max(distances))*100., "%")