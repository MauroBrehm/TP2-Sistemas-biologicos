import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

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
def viterbi(obs, estados, pi, a, b):
   '''Implementa el algoritmo de Viterbi para encontrar la secuencia de estados mas probable a partir
    de una secuencia de observaciones'''
   T = len(obs)
   estado_id = {s: i for i, s in enumerate(estados)}
   V = np.zeros((T, len(estados))) # Matriz de probabilidades V[t][s] de Viterbi
   psi = np.zeros((T, len(estados)), -1,dtype=int) # Matriz de estados anteriores psi[t][s]

   #Inicializacion
   for s in estados:
      i = estado_id[s]
      V[0][i] = pi[s] * b[s][obs[0]]
      #no hay estado anterior psi[0]

   #Recursion
   for t in range (1, T):
      for e in estados:
         i = estado_id[e]
         candidatos = [V[t-1][estado_id[er]] * a[er][e] for er in estados]
         psi[t][e] = np.argmax(candidatos) #Guardar el estado anterior s' que dio esa probabilidad máxima
         V[t][i] = max(candidatos) * b[e][obs[t]]
      
   #Al final (en t = T-1):
   best_final =np.argmax(V[T-1]) #Buscar el estado s con mayor V[T-1][s]
   prob = V[T-1][best_final] #probabilidad final prob = max(V[T-1][s])

   path = [0] * T
   path[T-1] = best_final
   for t in range(T-2, -1, -1):
      path[t] = psi[t+1][path[t+1]]

   path = [estados[i] for i in path] #Ruta óptima: path[s] correspondiente al estado más probable
   return path, prob # Devolver la ruta óptima (path) y su probabilidad (prob)

# ==================================================
#ACT 2 - Implemetacion Viterbi de decodificacion
# ==================================================
def reestimar_p(secuencia_obs, secuencia_estados, estados, simbolos):
   '''Reestima los parametros a partir de secuencia de observaciones y estados '''
   #Contar transiciones y emisiones
   cont_pi = defaultdict(float)
   cont_a = defaultdict(lambda: defaultdict(float))
   cont_b = defaultdict(lambda: defaultdict(float))
   cont_estados = defaultdict(float)

   for obs, path in zip(secuencia_obs, secuencia_estados):
      # pi: estados inicial
      cont_pi[path[0]] += 1

      for t in range(len(obs)):
         cont_b [path[t]][obs[t]] += 1 #contar emisiones
         cont_estados[path[t]] += 1 #contar estados

         if t < len(obs) - 1:
            cont_a[path[t]][path[t+1]] += 1 #contar transiciones
   #Reestimar parametros
   N = len(secuencia_obs)

   pi_n = {s: cont_pi[s] / N for s in estados}

   a_n = {}
   b_n = {}
   for s in estados:
      total_a = sum(cont_a[s][sp] for sp in estados) + len(estados)
      a_n[s] = {sp: (cont_a[s][sp] + 1)/ total_a for sp in estados}

      total_b = sum(cont_b[s][e] for e in simbolos) + len(simbolos)
      b_n [s] = {e: (cont_b[s][e] + 1) / total_b for e in simbolos}
   
   return pi_n, a_n, b_n

def decodificar_viterbi(secuencia_obs, estados, simbolos, max_iter = 100):
   '''Usa metodo iterativo de Viterbi '''

   #Inicializa parametros
   pi = {s: 1/len(estados) for s in estados}
   a ={s: {sp: 1/len(estados) for sp in estados} for s in estados}
   b = {s: {e: 1/len(simbolos) for e in simbolos} for s in estados}

   historial = []

   for iter in range(max_iter):
      paths = []
      log_prob_total = 0

      for obs in secuencia_obs:
         path, log_prob = viterbi(obs, estados, pi, a, b)
         paths.append(path)
         log_prob_total += log_prob
   historial.append(log_prob_total)

   #reestimar parametros
   pi_nuevo, a_nuevo, b_nuevo = reestimar_p(secuencia_obs, paths, estados, simbolos)

   #Verificamos la convergencia
   # No se como hacerlo :)

   return #se devuelve parametros finales y el historial de log_prob_total

# ==================================================
#ACT 3 - Secuencias de entrenamiento
# ==================================================
