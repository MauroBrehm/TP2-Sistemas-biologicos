import numpy as np
import matplotlib.pyplot as plt

#Paramentros del modelo real (lo que nos da el tp)

#Estados posibles del modelo oculto
estados=[ 'C', 'NC'] #C: codificante, NC: no codificante

#simbolos posibles de salida
simbolos=[ 'A', 'C', 'G', 'T']

#pi: probabilidad inicial de cada estado (con que probabilidad empieza la secuencia en cada estado)
#tp no dice uso 50/50
pi={'C': 0.5, 'NC': 0.5}

#a: matriz de transicion entre estados (probabilidad de pasar de un estado a otro)
a={'C': {'C': 0.995, 'NC': 0.005}, #desde C, 99.5% de las veces sigo en C, 0.5% paso a NC
   'NC': {'C': 0.010, 'NC': 0.990}} #desde NC, 1% paso a C, 99% sigo en NC

#b: matriz de emision (probabilidad de emitir cada simbolo dado el estado)
b={'C': {'A': 0.30, 'C': 0.20, 'G': 0.20, 'T': 0.30}, 
   'NC': {'A': 0.10, 'C': 0.40, 'G': 0.40, 'T': 0.10}} 

#==================================================
#ACT 1 - Algoritmo de Viterbi
#=====================================================
