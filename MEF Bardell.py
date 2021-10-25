# -*- coding: utf-8 -*-
"""
Created on Sun May 23 15:13:13 2021

@author: João Pedro Gomes
"""

#Importações
import math
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
from scipy import integrate
from Funçõesdeforma import func_forma

# Coordenadas da Barra
tamanho_barra = 337.7
tamanho_elemento = 337.7
tamanho_vetor =   math.ceil(tamanho_barra/tamanho_elemento)
coord = np.zeros((tamanho_vetor + 1,2))
vet_aux = 0
for i in range(tamanho_vetor):
    coord[i+1,1] = vet_aux + tamanho_elemento 
    vet_aux = vet_aux + tamanho_elemento
num_nos = coord.shape[0]

# Materiais (Modulo de elestacidade e Inercia)
mate = np.array([[20000,997]])

# Incidencia [nó inicial, nó final, material], um elemento por linha
num_mate = mate.shape[0]
if(num_mate == 1):
    num_ele = math.ceil(tamanho_barra/tamanho_elemento)
    inci = np.zeros((num_ele,3),dtype=int)
    for i in range(num_ele):
        k = i
        for j in range(2):
            inci[i][j] = k
            k = i + 1
    
else:
    inci = np.array([[0,1,0],[1,2,0]]) # Editar manualmente se for o caso de mais de um material ou seção
    num_ele = inci.shape[0]
    

# Plotagem da geometria

xn = np.transpose(coord[:,:1])
yn = np.transpose(coord[:,1:2])
plt.scatter(xn,yn,color="red")
for i in range(num_ele):
    xe = np.array([coord[inci[i,0],0],coord[inci[i,1],0]])
    ye = np.array([coord[inci[i,0],1],coord[inci[i,1],1]])
    plt.plot(xe,ye,color="orange")
    
plt.title("Malha da Barra")
plt.show()

# MEF

# MEF Matrizes Globais

num_eqs = 10
k_ele = np.zeros((num_eqs,num_eqs))
kg_ele = np.zeros((num_eqs,num_eqs))


# MEF - Matriz de Rigidez

keg = np.zeros((2*num_nos + (num_eqs-4)*num_ele,2*num_nos + (num_eqs-4)*num_ele))
for i in range(num_ele):
    no_inic = inci[i,0]
    no_final = inci[i,1]
    xni = coord[no_inic,0]
    yni = coord[no_inic,1]
    xnf = coord[no_final,0]
    ynf = coord[no_final,1]
    le = ((xnf-xni)**2 + (ynf-yni)**2)**(1/2)
    E = mate[inci[i,2],0]
    I = mate[inci[i,2],1]
    ff = func_forma(le)
    dif2 = ff.dif2
    for v in range(num_eqs):
        for b in range(num_eqs):
            produto = lambda x: dif2[v](x)*dif2[b](x)
            k_ele[v][b] = ((E*I)*8/le**3)*float(integrate.fixed_quad(produto,-1,1, n =100)[0])
    posi = [2*no_inic,2*no_inic+1,2*no_final,2*no_final+1,(2*num_ele+2)+(6*i),(2*num_ele+3)+(6*i),(2*num_ele+4)+(6*i),(2*num_ele+5)+(6*i),(2*num_ele+6)+(6*i),(2*num_ele+7)+(6*i)]
    for j in range(num_eqs):
        for k in range(num_eqs):
            keg[posi[j],posi[k]] += k_ele[j,k]

# Matriz Geometrica
kg = np.zeros((2*num_nos + (num_eqs-4)*num_ele,2*num_nos + (num_eqs-4)*num_ele))
for i in range(num_ele):
    no_inic = inci[i,0]
    no_final = inci[i,1]
    xni = coord[no_inic,0]
    yni = coord[no_inic,1]
    xnf = coord[no_final,0]
    ynf = coord[no_final,1]
    le = ((xnf-xni)**2 + (ynf-yni)**2)**(1/2)
    ff = func_forma(le)
    dif1 = ff.dif1
    for v in range(num_eqs):
        for b in range(num_eqs):
           produto = lambda x: dif1[v](x)*dif1[b](x)
           kg_ele[v][b] = (2/le)*integrate.fixed_quad(produto,-1,1, n = 100)[0]
    posi = [2*no_inic,2*no_inic+1,2*no_final,2*no_final+1,(2*num_ele+2)+(6*i),(2*num_ele+3)+(6*i),(2*num_ele+4)+(6*i),(2*num_ele+5)+(6*i),(2*num_ele+6)+(6*i),(2*num_ele+7)+(6*i)]
    for j in range(num_eqs):
        for k in range(num_eqs):
            kg[posi[j],posi[k]] += kg_ele[j,k]


# Condições de Contorno Essenciais
num_cond = 3
ordem = [3,4,5,6,7,8,9,0,1,2]

# Reordenação das Matrizes
kg_ord = np.zeros((2*num_nos + (num_eqs-4)*num_ele,2*num_nos + (num_eqs-4)*num_ele))
keg_ord = np.zeros((2*num_nos + (num_eqs-4)*num_ele,2*num_nos + (num_eqs-4)*num_ele))
for i in range(2*num_nos + (num_eqs-4)*num_ele):
    for j in range(2*num_nos + (num_eqs-4)*num_ele):
        kg_ord[i,j] = kg[ordem[i],ordem[j]]
        keg_ord[i,j] = keg[ordem[i],ordem[j]]

# Solução
dimkll = 2*num_nos + (num_eqs-4)*num_ele - num_cond
P = linalg.eigh(keg_ord[:dimkll,:dimkll],kg_ord[:dimkll,:dimkll])
print("P = {}".format(P[0][0]))

#Plotagem do grafico

#Reordenar o autovalor na sequencia dos graus de liberadade
autov = np.zeros((2*num_nos + 6*num_ele,1))
for i in range(dimkll):
    autov[ordem[i],0] = P[1][i,0]

# Pós-processamento
for i in range(num_ele):
    no_inic = inci[i,0]
    no_final = inci[i,1]
    xni = coord[no_inic,0]
    yni = coord[no_inic,1]
    xnf = coord[no_final,0]
    ynf = coord[no_final,1]
    le = ((xnf-xni)**2 + (ynf-yni)**2)**(1/2)
    E = mate[inci[i,2],0]
    I = mate[inci[i,2],1]
    posi = [2*no_inic,2*no_inic+1,2*no_final,2*no_final+1,(2*num_ele+2)+(6*i),(2*num_ele+3)+(6*i),(2*num_ele+4)+(6*i),(2*num_ele+5)+(6*i),(2*num_ele+6)+(6*i),(2*num_ele+7)+(6*i)]
    valid_plot = np.zeros((10))
    def n1(x):
        ff = func_forma(le)
        hermite = ff.hermite
        for j in range(num_eqs):
            valid_plot[j] = 1
        return (hermite[0](x)*autov[posi[0],0]*valid_plot[0]) + (hermite[1](x)*autov[posi[1],0]*valid_plot[1]) + (hermite[2](x)*autov[posi[2],0]*valid_plot[2]) + (hermite[3](x)*autov[posi[3],0]*valid_plot[3]) + (hermite[4](x)*autov[posi[4],0]*valid_plot[4]) + (hermite[5](x)*autov[posi[5],0]*valid_plot[5]) + (hermite[6](x)*autov[posi[6],0]*valid_plot[6]) + (hermite[7](x)*autov[posi[7],0]*valid_plot[7]) + (hermite[8](x)*autov[posi[8],0]*valid_plot[8]) + (hermite[9](x)*autov[posi[9],0]*valid_plot[9])
        #return ((((1/2)-(3/4)*(((2*x)/le)-1)+(1/4)*(((2*x)/le)-1)**3)*autov[posi[0],0]) + ((le/8)*(1-(((2*x)/le)-1)-(((2*x)/le)-1)**2+(((2*x)/le)-1)**3)*autov[posi[1],0]) + (((1/2)+(3/4)*(((2*x)/le)-1)-(1/4)*(((2*x)/le)-1)**3)*autov[posi[2],0]) + ((le/8)*(-1-(((2*x)/le)-1)+(((2*x)/le)-1)**2+(((2*x)/le)-1)**3)*autov[posi[3],0]))
    x = np.linspace(0,le,num=100)
    plt.plot(n1(x),x+(i*le),color="blue")
xn = np.transpose(coord[:,:1])
yn = np.transpose(coord[:,1:2])
plt.scatter(xn,yn,color="red")
for i in range(num_ele):
    xe = np.array([coord[inci[i,0],0],coord[inci[i,1],0]])
    ye = np.array([coord[inci[i,0],1],coord[inci[i,1],1]])
    plt.plot(xe,ye,color="orange")
plt.xticks(range(-5, 5))
plt.title("Modo de Flambagem")
plt.show()