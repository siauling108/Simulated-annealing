# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 09:53:20 2020

@author: marco
"""
import matplotlib.pyplot as plt
import configparser
#import sys
#from sys import argv
import functions


config = configparser.ConfigParser()
config.read('input.txt')
#config.read(sys.argv[1])

N = config.getint('parameters', 'N') 
iterations = config.getint('parameters', 'iterations')
nstep = config.getint('parameters', 'nstep')
T_min = config.getfloat('parameters', 'T_min')
T = config.getfloat('parameters', 'T')
alpha = config.getfloat('parameters', 'alpha')

#Tem = []
#accelist = []
city = functions.travel(N)
citynum = list(range(N))
start=functions.length(N, city, citynum)
print('Initial distance:',start)
distances=[start]
functions.plpath(N, city, citynum)

'''
OPTIONS (RAW)
Moves: swap, block reverse, prune and graft
Criteria: Metropolis, distance
'''

for ii in range(iterations):
    #citynum, accelist, Tem=functions.anneal_BRev_distance(N, alpha, T, city, T_min, nstep, citynum)
    #citynum, accelist, Tem=functions.anneal_BRev_Metropolis(N, alpha, T, city, T_min, nstep, citynum) 
    #citynum, accelist, Tem=functions.anneal_swap_Metropolis(N, alpha, T, city, T_min, nstep, citynum) #BOH
    citynum, accelist, Tem=functions.anneal_swap_distance(N, alpha, T, city, T_min, nstep, citynum) 
    #citynum, accelist, Tem=functions.anneal_PG_Metropolis(N, alpha, T, city, T_min, nstep, citynum)
    #citynum, accelist, Tem=functions.anneal_PG_distance(N, alpha, T, city, T_min, nstep, citynum)

    if ii==0:
        functions.acceptance_plot(accelist, Tem)
    distances.append(functions.length(N, city, citynum))
    print('Current distance:',distances[ii+1],', Iteration:',ii+1)
    functions.plpath(N, city, citynum)
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
print("Percentage decrease of the distance", (distances[0]-distances[-1])/distances[0]*100., "%")
#print("Percentage decrease of the distance", (min(distances)/max(distances))*100., "%")