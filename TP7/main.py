import numpy as np
import matplotlib.pyplot as plt

#Parametros del modelo
filas=50
cols=50 #estos valores los puse yo random no dice que valore tiene que tener
dt=1.0 #paso de tiempo en milisegundos

#umbrales de exitacion
UR= 90.0 #umbral para pasar de reposo a excitacion
UP= 120.0 #umbral para reexitacion desde PRR

#Conductividad entre celulas vecinas
G = 1.0

#Umbrales para cambios de estado
v_reposo = -90.0
v_pico= 30.0
v_prr = -0.001

#tasas de cambio 
t_existacion =60.0 #mV/ms --> velocidad de subida en estado E
t_repolariz=-1.0 #mV/ms --> velocidad de bajada en estado PRR
t_pra_exp =0.04 #fraccion por ms --> decaimiento exponencial en PRA (4% por ms)

#Marcapasos --> frecuencia cardíaca de 60 latidos/minuto = 1 latido cada 1000 ms
periodo_marcapasos = 1000.0 #ms

#Duracion total de la simulacion
duracion = 5000 #ms (5 segundos --> 5 latidos)
n_pasos = int(duracion / dt)

#------------------------------------------------------------------------------
#Defincion de estados 
reposo = 0 #R --> celula en reposo
excitado = 1 #E --> celula excitada
prr = 3 #periodo refractario relativo --> puede re-excitarse con estimulo fuerte
pra = 2 #periodo refractario absoluto --> no puede re-excitarse
#------------------------------------------------------------------------------
#Inicializacion de variables de estado

#Matriz de estados --> todas las celulas empiezan en reposo
estado = np.full((filas, cols), reposo, dtype=int) 

#matriza de potenciales de membrana --> todas las celulas empiezan en potencial de reposo
potencial = np.full((filas, cols), v_reposo, dtype=float) 

#registros de ECG a lo largo del tiempo
ECG =[]

