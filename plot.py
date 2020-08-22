# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 10:51:34 2020

@author: marco
"""
import numpy as np
import matplotlib.pyplot as plt
import configparser

config = configparser.ConfigParser()
config.read('inpus.txt')

iterations = config.getint('parameters', 'iterations')
file1 = config.get('files','distances')
file2 = config.get('files','path')
file3 = config.get('files','Tem')
file4 = config.get('files','tot_accept')

plot_a = config.get('files','path_pl')
plot_b = config.get('files','acc_pl')
plot_c = config.get('files','optimization_pl')


def dist_optimization():
    '''
    Plots the travel length as a function of the number of iterations.
    '''
    distances = np.load(file1)
    fig = plt.figure()
    plt.title("Current distance vs. # of iterations")
    plt.xlabel("Iteration")
    plt.ylabel("Current distance (arb. units)")
    plt.plot(distances)
    plt.grid()
    plt.show()   
    
    fig.savefig(plot_c)
    
    
def path_plot():
    '''
    Plots the followed path.
    '''
    
    path = np.load(file2)
    
    fig = plt.figure()
    plt.title("Followed path")
    plt.plot (path[0],path[1],'o-')
    plt.xlabel("x coordinate (arb. units)")
    plt.ylabel("y coordinate (arb. units)")    
    plt.grid()
    plt.show()
    
    fig.savefig(plot_a)

 
    
def acceptance_plot(iterations):
    '''
    Plots the acceptance rate as a function of the temperature for each iteration.
    '''
    Tem = np.load(file3)
    tot_acceptance = np.load(file4)
    fig = plt.figure()    
    plt.title("Moves acceptance rate vs. temperature")
    plt.xlabel("Temperature (arb. units)")
    plt.ylabel("Moves acceptance rate (%)") 
    lab=[]
    for i in range(iterations): 
        plt.plot(Tem,tot_acceptance[i])
        lab.append(i+1) 
    plt.legend(lab)
    plt.grid()
    plt.show()
    
    fig.savefig(plot_b)


dist_optimization()
path_plot()
acceptance_plot(iterations)