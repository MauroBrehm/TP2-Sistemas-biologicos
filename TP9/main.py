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
pi={'C': 0.6, 'NC': 0.4}

#a: matriz de transicion entre estados (probabilidad de pasar de un estado a otro)
a={'C': {'C': 0.995, 'NC': 0.005}, #desde C, 99.5% de las veces sigo en C, 0.5% paso a NC
   'NC': {'C': 0.010, 'NC': 0.990}} #desde NC, 1% paso a C, 99% sigo en NC

#b: matriz de emision (probabilidad de emitir cada simbolo dado el estado)
b={'C': {'A': 0.30, 'C': 0.20, 'G': 0.20, 'T': 0.30}, 
   'NC': {'A': 0.10, 'C': 0.40, 'G': 0.40, 'T': 0.10}} 

# =============================================================================
# LECTURA DEL ARCHIVO CSV
# =============================================================================
# Columna 'estados': 1 = C (codificante), 2 = NC (no codificante)
# Columna 'salidas': 1=A, 2=T, 3=C, 4=G  
 
NUM_A_SIMBOLO = {1: 'A', 2: 'T', 3: 'C', 4: 'G'}
NUM_A_ESTADO  = {1: 'C', 2: 'NC'}
 
def cargar_csv(ruta):
    """
    Lee el CSV línea por línea y devuelve:
      obs_dec          : lista de bases de la secuencia de decodificación
      est_dec_real     : lista de estados reales de la secuencia de decodificación
      secs_ap          : lista de 4 listas de bases (secuencias de aprendizaje)
      secs_ap_estados  : lista de 4 listas de estados reales (aprendizaje)
    """
    # Usamos un dict para ir juntando cada secuencia por su nombre
    obs_por_indice    = defaultdict(list)   # bases
    est_por_indice    = defaultdict(list)   # estados reales
    tipo_por_indice   = {}                  # qué tipo es cada índice
 
    with open(ruta, 'r') as f:
        lineas = f.readlines()

    for linea in lineas[1:]:
        partes = linea.strip().split(',')
        tipo    = partes[0]   
        indice  = partes[1]   
        estado  = int(partes[2])   
        salida  = int(partes[3])   
        
        obs_por_indice[indice].append(NUM_A_SIMBOLO[salida])
        est_por_indice[indice].append(NUM_A_ESTADO[estado])
        tipo_por_indice[indice] = tipo
 
    # --- Separar secuencia de decodificación ---
    obs_dec      = obs_por_indice['sec_dec_1']
    est_dec_real = est_por_indice['sec_dec_1']
 
    # --- Separar secuencias de aprendizaje (en orden: ap_1, ap_2, ap_3, ap_4)----
    indices_ap = sorted(k for k, v in tipo_por_indice.items()
                        if v == 'sec_aprendizaje')
 
    secs_ap         = [obs_por_indice[idx] for idx in indices_ap]
    secs_ap_estados = [est_por_indice[idx] for idx in indices_ap]
 
    return obs_dec, est_dec_real, secs_ap, secs_ap_estados

#==================================================
#ACT 1 - Algoritmo de Viterbi
#=====================================================
def viterbi(obs, estados, pi, a, b):
   '''Implementa el algoritmo de Viterbi para encontrar la secuencia de estados mas probable a partir
    de una secuencia de observaciones'''
   T = len(obs)
   N = len(estados)

   #Mapa de nombre de estado -> indice numerico para indexar arrayd
  #Ej: {'C': 0, 'NC': 1}
   estado_id = {s: i for i, s in enumerate(estados)}

   #V[t][i]= log-probabilidad del camino mas probable hasta el tiempo t terminando en el estado i
   V = np.full((T, N), -np.inf) # Matriz de probabilidades V[t][s] de Viterbi
   psi = np.full((T, N), -1, dtype=int) # Matriz de estados anteriores psi[t][s]

   #Inicializacion (t=0)
   #V[0][i]= log(pi[s]) + log(b[s][obs[0]]) para cada estado s
   for s in estados:
      i = estado_id[s]
      p_pi=pi[s]
      p_b=b[s][obs[0]]
      #solo calculamos si ambas probabilidades son positivas
      if p_pi > 0 and p_b > 0:
         V[0][i] = np.log(p_pi) + np.log(p_b) 
         #Si alguna de las probabilidades es 0, dejamos V[0][i] como -inf (log(0)=-inf) para indicar que esa ruta no es posible
      
      #no hay estado anterior psi[0]

   #Recursion (t=1 a T-1)
   #V[t][i] = max_{j} (V[t-1][j] + log(a[j][i]) + log(b[i][obs[t]])) para cada estado i
   for t in range (1, T):
      for e in estados:
         i = estado_id[e] #indice del estado actual e
         p_b = b[e][obs[t]]
         if p_b == 0:
            continue #Si la probabilidad de emisión es 0, esa ruta no es posible, dejamos V[t][i] como -inf
         log_p_b = np.log(p_b)
         #Calculamos la probabilidad de llegar al estado e desde cada estado anterior posible y nos quedamos con la mayor
         mejor_prob = -np.inf
         mejor_prev = -1

         for e_prev in estados:
            j = estado_id[e_prev] #indice del estado ANTERIOR e_prev
            p_a = a[e_prev][e]
            if p_a == 0:
               continue #Si la probabilidad de transición es 0, esa ruta no es posible, seguimos con el siguiente estado anterior
            candidato = V[t-1][j] + np.log(p_a)
            if candidato > mejor_prob:
               mejor_prob = candidato
               mejor_prev = j
         V[t][i] = mejor_prob + log_p_b #Probabilidad total de la mejor ruta hasta el estado e en el tiempo t
         psi[t][i] = mejor_prev #Guardamos el estado anterior que nos dio la mejor probabilidad
      
   #Al final (en t = T-1):
   best_final =int(np.argmax(V[T-1])) #Buscar el estado s con mayor V[T-1][s]
   prob = V[T-1][best_final] #probabilidad final prob = max(V[T-1][s])

   path = [0] * T
   path[T-1] = best_final
   for t in range(T-2, -1, -1):
      #el estado en t es el predecesor del estado en t+1
      path[t] = psi[t+1][path[t+1]]

   path = [estados[i] for i in path] #Ruta óptima: path[s] correspondiente al estado más probable
   return path, prob # Devolver la ruta óptima (path) y su probabilidad (prob)

# ==================================================
#ACT 2 - Implemetacion Viterbi de decodificacion y ACT 3 - secuencia de entrenamiento
# ==================================================
def reestimar_p(secuencia_obs, secuencia_estados, estados, simbolos):
   '''Reestima los parametros a partir de secuencia de observaciones y estados '''
   #Contar transiciones y emisiones
   cont_pi = defaultdict(float) #Contador de estados iniciales
   cont_a = defaultdict(lambda: defaultdict(float)) #Contador de transiciones entre estados
   cont_b = defaultdict(lambda: defaultdict(float)) #Contador de emisiones
   N_secuencias = len(secuencia_obs) #Cantidad de secuencias de entrenamiento

   for obs, path in zip(secuencia_obs, secuencia_estados):
      # Estado inicial de esta secuencia
      cont_pi[path[0]] += 1

      for t in range(len(obs)):
         #contar la emision, en el estado path[t] se emitio obs[t]
         cont_b [path[t]][obs[t]] += 1 #contar emisiones

         #contar la transicion (si no es el ultimo paso)
         if t < len(obs) - 1:
            cont_a[path[t]][path[t+1]] += 1 #contar transiciones
   #Reestimar parametros
   #pi[s]= veces que la secuencia empieza en s / total de secuencias
   pi_n = {s: cont_pi[s] / N_secuencias for s in estados}

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
   paths_finales = None

   for iter in range(max_iter):
      paths = []
      log_prob_total = 0

      for obs in secuencia_obs:
         path, log_prob = viterbi(obs, estados, pi, a, b)
         paths.append(path)
         log_prob_total += log_prob
      historial.append(log_prob_total)
      paths_finales = paths

   #reestimar parametros
   pi_nuevo, a_nuevo, b_nuevo = reestimar_p(secuencia_obs, paths, estados, simbolos)

def entrenar_viterbi(secuencias_obs, estados, simbolos, max_iter=100, tolerancia=1e-4, tol=None): #tolerancia es la tolerancia para la convergencia, 
   #cuahndo ya casi no cambia  la log-probabilidad total entre iteraciones, consideramos que el modelo ha convergido y cortamos el entrenamiento
    """
    Entrena el HMM iterando Viterbi + re-estimacion hasta convergencia
    tolerancia: si la log-prob total cambia menos que este valor entre dos iteraciones consecutivas, consideramos que el modelo convergio y cortamos
    """
    if tol is not None:
        tolerancia = tol
    pi_cur = {s: 1/len(estados) for s in estados}
    a_cur  = {s: {sp: 1/len(estados) for sp in estados} for s in estados}
    b_cur  = {s: {e:  1/len(simbolos) for e  in simbolos} for s in estados}
 
    historial     = []
    paths_finales = None
 
    for it in range(max_iter):
        paths          = []
        log_prob_total = 0.0
 
        for obs in secuencias_obs:
            path, lp = viterbi(obs, estados, pi_cur, a_cur, b_cur)
            paths.append(path)
            log_prob_total += lp
 
        historial.append(log_prob_total)
        paths_finales = paths
 
        pi_nuevo, a_nuevo, b_nuevo = reestimar_p(
            secuencias_obs, paths, estados, simbolos
        )
 
        if it > 0 and abs(historial[-1] - historial[-2]) < tolerancia:
            print(f"  Convergió en iteración {it+1}  "
                  f"(log-prob = {abs(historial[-1]-historial[-2]):.2e})")
            pi_cur, a_cur, b_cur = pi_nuevo, a_nuevo, b_nuevo
            break
 
        pi_cur, a_cur, b_cur = pi_nuevo, a_nuevo, b_nuevo
    else:
        print(f"  Alcanzó {max_iter} iteraciones sin converger.")
 
    return pi_cur, a_cur, b_cur, historial, paths_finales



# =============================================================================
#Grafica
# =============================================================================
def graficar_convergencia(historial):
    plt.figure(figsize=(7, 4))
    plt.plot(range(1, len(historial)+1), historial,
             marker='o', markersize=4, linewidth=1.5, color='steelblue')
    plt.xlabel("Iteración")
    plt.ylabel("Log-probabilidad total")
    plt.title("Convergencia del entrenamiento (Viterbi iterativo)")
    plt.tight_layout()
    plt.show()

# =============================================================================
#archivo y funciones
# =============================================================================
import os

# __file__ es la ruta del script actual
# dirname le saca el nombre del archivo y se queda con la carpeta
# join une la carpeta con el nombre del csv dentro de modules
RUTA_CSV = os.path.join(os.path.dirname(__file__), 'modules', 'secuencias.csv')

obs_dec, est_dec_real, secs_ap, secs_ap_estados = cargar_csv(RUTA_CSV)
print(f"\nDatos cargados:")
print(f"  Secuencia de decodificación: {len(obs_dec)} bases")
print(f"  Secuencias de aprendizaje:   {len(secs_ap)} × {len(secs_ap[0])} bases")

print("ACT 1 — Viterbi con parámetros reales (sec_decodificacion)")
print("-"*58)
path_dec, prob_dec = viterbi(obs_dec, estados, pi, a, b)

# --- Actividades 2 y 3 ---
print("\n" + "-"*58)
print("ACT 2 y 3 — Entrenamiento con sec_ap_1 a sec_ap_4")
print("-"*58)
 
pi_est, a_est, b_est, historial, paths_ap = entrenar_viterbi(
      secs_ap, estados, simbolos, max_iter=100, tolerancia=1e-4
    )
graficar_convergencia(historial)

#Actividad 4: Comparamos los parametros del sistema contra los parametros estimados 
def comparar_parametros(pi_real, a_real, b_real,pi_est, a_est, b_est,estados, simbolos):
   print("\n" + "="*60)
   print("ACT 4 — COMPARACIÓN DE PARÁMETROS")
   print("="*60)

    
   # Probabilidades iniciales
    
   print("\nProbabilidades iniciales (π)")
   print("-"*40)

   for s in estados:
      real = pi_real[s]
      est  = pi_est[s]
      error = abs(real - est)

      print(f"{s:>3} | Real = {real:.4f} | "
            f"Estimado = {est:.4f} | "
            f"Error = {error:.4f}")
   # Matriz de transición
   print("\nMatriz de transición (A)")
   print("-"*40)
   for s in estados:
        for sp in estados:
            real = a_real[s][sp]
            est  = a_est[s][sp]
            error = abs(real - est)

            print(f"{s}->{sp} | Real = {real:.4f} | "
                  f"Estimado = {est:.4f} | "
                  f"Error = {error:.4f}")
   # Matriz de emisión
   print("\nMatriz de emisión (B)")
   print("-"*40)

   for s in estados:
        for e in simbolos:
            real = b_real[s][e]
            est  = b_est[s][e]
            error = abs(real - est)

            print(f"{s} emite {e} | "
                  f"Real = {real:.4f} | "
                  f"Estimado = {est:.4f} | "
                  f"Error = {error:.4f}")

comparar_parametros(pi, a, b, pi_est, a_est, b_est, estados, simbolos)
#Graficamos la comparacion entre lo que nos dio y lo que deberia darnos
def graficar_comparacion(pi_real, a_real, b_real, pi_est, a_est, b_est, estados, simbolos):
    # Gráfica de comparación de parámetros
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # Probabilidades iniciales
    axes[0].bar(estados, [pi_real[s] for s in estados], alpha=0.5, label='Real')
    axes[0].bar(estados, [pi_est[s] for s in estados], alpha=0.5, label='Estimado')
    axes[0].set_title("Probabilidades iniciales (π)")
    axes[0].legend()

    # Matriz de transición
    a_real_vals = [[a_real[s][sp] for sp in estados] for s in estados]
    a_est_vals  = [[a_est[s][sp] for sp in estados] for s in estados]

    im1 = axes[1].imshow(a_real_vals, cmap='Blues', vmin=0, vmax=1)
    im2 = axes[1].imshow(a_est_vals, cmap='Reds', vmin=0, vmax=1, alpha=0.5)
    axes[1].set_xticks(range(len(estados)))
    axes[1].set_yticks(range(len(estados)))
    axes[1].set_xticklabels(estados)
    axes[1].set_yticklabels(estados)
    axes[1].set_title("Matriz de transición (A)")
    fig.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    # Matriz de emisión
    b_real_vals = [[b_real[s][e] for e in simbolos] for s in estados]
    b_est_vals  = [[b_est[s][e] for e in simbolos] for s in estados]

    im3 = axes[2].imshow(b_real_vals, cmap='Blues', vmin=0, vmax=1)
    im4 = axes[2].imshow(b_est_vals, cmap='Reds', vmin=0, vmax=1, alpha=0.5)
    axes[2].set_xticks(range(len(simbolos)))
    axes[2].set_yticks(range(len(estados)))
    axes[2].set_xticklabels(simbolos)
    axes[2].set_yticklabels(estados)
    axes[2].set_title("Matriz de emisión (B)")
    fig.colorbar(im3, ax=axes[2], fraction=0.046, pad=0.04)

graficar_comparacion(pi, a, b, pi_est, a_est, b_est, estados, simbolos)

#Probando cosas
print("\n¡A ver si anda")

for i, path in enumerate(paths_ap):
    nC = path.count('C')
    nNC = path.count('NC')

    print(f"Secuencia {i+1}")
    print(f"C  = {nC}")
    print(f"NC = {nNC}")
print("Fin")
pi_real_obs, a_real_obs, b_real_obs = reestimar_p(
    secs_ap,
    secs_ap_estados,
    estados,
    simbolos
)
comparar_parametros(pi, a, b, pi_real_obs, a_real_obs, b_real_obs, estados, simbolos)

