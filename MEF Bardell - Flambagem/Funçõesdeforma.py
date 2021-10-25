# -*- coding: utf-8 -*-
"""
Created on Sun May 23 15:13:13 2021

@author: João Pedro Gomes
"""
import numpy as np

class func_forma:
    def __init__(self, le):
        # Funções de forma

        h1 = lambda x: 1/2 - (3/4)*(((2*x)/le)-1) + (1/4)*(((2*x)/le)-1)**3
        h2 = lambda x: (le/8)*(1- (((2*x)/le)-1) - (((2*x)/le)-1)**2+ (((2*x)/le)-1)**3)
        h3 = lambda x: 1/2 + (3/4)*(((2*x)/le)-1) - (1/4)*(((2*x)/le)-1)**3
        h4 = lambda x: (le/8)*(-1 - (((2*x)/le)-1) + (((2*x)/le)-1)**2  + (((2*x)/le)-1)**3)
        h5 = lambda x: 1/8 - (1/4)*((((2*x)/le)-1)**2) + (1/8)*(((2*x)/le)-1)**4
        h6 = lambda x: (1/8)*(((2*x)/le)-1) - (1/4)*((((2*x)/le)-1)**3) + (1/8)*(((2*x)/le)-1)**5
        h7 = lambda x: -1/48 + (3/16)*((((2*x)/le)-1)**2) - (5/16)*((((2*x)/le)-1)**4) + (7/48)*(((2*x)/le)-1)**6
        h8 = lambda x: -(1/16)*(((2*x)/le)-1) + (5/16)*((((2*x)/le)-1)**3) - (7/16)*((((2*x)/le)-1)**5) + (3/16)*(((2*x)/le)-1)**7
        h9 = lambda x: 1/128 - (5/32)*((((2*x)/le)-1)**2) + (35/64)*((((2*x)/le)-1)**4) - (21/32)*((((2*x)/le)-1)**6) + (33/128)*(((2*x)/le)-1)**8
        h10 = lambda x: (5/128)*(((2*x)/le)-1) - (35/96)*((((2*x)/le)-1)**3) + (63/64)*((((2*x)/le)-1)**5) - (33/32)*((((2*x)/le)-1)**7) + (143/384)*(((2*x)/le)-1)**9
        hermite = np.array((h1,h2,h3,h4,h5,h6,h7,h8,h9,h10))
        self.hermite = hermite
        
        # Primeira derivada
        g1 = lambda x: -3/4 + (3/4)*x**2
        g2 = lambda x: (le/8)*(3*x**2 - 2*x - 1)
        g3 = lambda x: 3/4 -(3/4)*x**2
        g4 = lambda x: (le/8)*(3*x**2 + 2*x - 1)
        g5 = lambda x: -(1/2)*x + (1/2)*(x**3)
        g6 = lambda x: 1/8 - (3/4)*(x**2) + (5/8)*x**4
        g7 = lambda x: (3/8)*x - (5/4)*(x**3) + (7/8)*x**5
        g8 = lambda x: -1/16 + (15/16)*(x**2) - (35/16)*(x**4) + (21/16)*x**6
        g9 = lambda x: -(5/16)*x + (35/16)*(x**3) - (63/16)*(x**5) + (33/16)*x**7
        g10 = lambda x: 5/128 - (35/32)*(x**2) + (315/64)*(x**4) - (231/32)*(x**6) + (429/128)*x**8
        dif1 = np.array((g1,g2,g3,g4,g5,g6,g7,g8,g9,g10))
        self.dif1 = dif1
         
         # Segunda Derivada   
        f1 = lambda x: (3/2)*x 
        f2=  lambda x: (le/8)*(6*x-2)
        f3 = lambda x: -(3/2)*x
        f4 = lambda x: (le/8)*(6*x+2)
        f5 = lambda x: -1/2 + (3/2)*x**2
        f6 = lambda x: -(3/2)*x + (5/2)*x**3
        f7 = lambda x: 3/8 - (15/4)*(x**2) + (35/8)*x**4
        f8 = lambda x: (15/8)*x - (35/4)*(x**3) + (63/8)*x**5
        f9 = lambda x: -5/16 + (105/16)*(x**2) - (315/16)*(x**4) + (231/16)*x**6
        f10= lambda x: -(35/16)*x + (315/16)*(x**3) - (693/16)*(x**5) + (429/16)*x**7
        dif2 = np.array((f1,f2,f3,f4,f5,f6,f7,f8,f9,f10))
        self.dif2 = dif2
        
        pass
    
