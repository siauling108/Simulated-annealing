# -*- coding: utf-8 -*-
"""
Created on Wed Aug 19 09:57:11 2020

@author: marco
"""
import numpy as np
from hypothesis import strategies as st
from hypothesis import given, assume
import hypothesis.extra.numpy as xn
import functions
import configparser


config = configparser.ConfigParser()
config.read('inputs.txt')

N = config.getint('parameters', 'N') 
iterations = config.getint('parameters', 'iterations')
nstep = config.getint('parameters', 'nstep')
T_min = config.getfloat('parameters', 'T_min')
T = config.getfloat('parameters', 'T')
alpha = config.getfloat('parameters', 'alpha')


@given(M=st.integers(2,N))
def test_leng_positive(M): 
    '''
    Tests if the travel lenght is actually positive.
    '''
    city = functions.travel(M)
    citynum = list(range(M))
    assume(len(city[:]==len(citynum)))
    distance = functions.length(M, city, citynum)
    assert distance > 0
    
    
@given(x=st.integers(0,N-1), y=st.integers(0,N-1), citynum=st.lists(st.integers(1, N-1),min_size=2,max_size=N,unique=True))   
def test_breverse(x, y, citynum):  
    '''
    Checks if the block reverse move works correctly.
    '''
    assume(x != y)
    citynum1=functions.breverse(x,y, citynum)
    citynum2=functions.breverse(y,x, citynum1)
    assert citynum1 == citynum2

    
@given(x=st.integers(0,N-1), y=st.integers(0,N-1), citynum=st.lists(st.integers(1, N-1),min_size=2,max_size=N,unique=True))   
def test_swap(x, y, citynum):  
    '''
    Checks if the swap move works correctly.
    '''
    assume(x != y)
    citynum1=functions.swap(x,y, citynum)
    citynum2=functions.swap(y,x, citynum1)
    assert citynum1 == citynum2


