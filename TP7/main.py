import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap

from modules.graficador import graficar_resultados


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
v_prr = 0.001

# Tasas de cambio
t_exc = 60.0   # mV/ms --> velocidad de subida en estado excitado

t_rep = 1.0    # mV/ms --> velocidad de bajada en PRR

t_pra_exp = 0.04  # fracción de decaimiento en PRA por ms

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


def potencial_total_vertical(matriz):
    """Calcula el potencial total vertical como suma de los dipolos verticales."""
    dipolos = matriz[:-1, :] - matriz[1:, :]
    return np.sum(dipolos)

#registros de ECG a lo largo del tiempo
ECG =[]
potencial_celula=[]
#------------------------------------------------------------------------+
#Graficaion de cuadrilla
# Colores para cada estado
# 0 = reposo
# 1 = excitado
# 2 = PRA
# 3 = PRR

cmap = ListedColormap([
    'blue',    # reposo
    'red',     # excitado
    'yellow',  # PRA
    'green'    # PRR
])

plt.ion()

fig, ax = plt.subplots(figsize=(6,6))
legendas = [
    mpatches.Patch(color='blue', label='Reposo'),
    mpatches.Patch(color='red', label='Excitado'),
    mpatches.Patch(color='yellow', label='PRA'),
    mpatches.Patch(color='green', label='PRR')
]

ax.legend(handles=legendas,
          bbox_to_anchor=(1.05,1),
          loc='upper left')
#------------------------------------------------------------------------------
for paso in range(n_pasos):
    if paso % int(periodo_marcapasos/dt) == 0:
        if estado[0, cols//2] == reposo: #solo se excita si esta en reposo, sino no se excita y se espera al siguiente periodo del marcapasos
        #se excita una sola celula en el centro de la matriz
            potencial[0, cols//2] = v_reposo
            estado[0, cols//2] = excitado

    nuevo_potencial = potencial.copy()
    nuevo_estado = estado.copy()

    for i in range(filas):
        for j in range(cols):
            V = potencial[i,j]
            e = estado[i,j]

            #Tema de los vecinos
            vecinos = []
            if i > 0:
                vecinos.append(potencial[i-1,j]) #vecino de arriba
            if i < filas-1:
                vecinos.append(potencial[i+1,j]) #vecino de abajo
            #lateral izquierdo periodico --> el vecino de la izquierda del borde se conecta con el de la derecha
            vecinos.append(potencial[i,(j-1)%cols])
            #lateral derecho periodico
            vecinos.append(potencial[i,(j+1)%cols]) #vecino de derecha

            I = G*sum(v - V for v in vecinos) #Corriente total recibida

            # --- Reposo -> Excitado ---
            if e == reposo and I >= UR:
                nuevo_estado[i,j] = excitado
                nuevo_potencial[i,j] = v_reposo

            elif e == excitado:     #Entra sodio por lo que el potencial sube rapido(ocurre despolarizacion)
                nuevo_potencial[i,j] = min(v_pico, V + t_exc * dt)
                # --- Excitado -> PRA ---
                if nuevo_potencial[i,j] >= v_pico: #cuando alcanza el pico, pasa a PRA
                    nuevo_estado[i,j] = pra
                    nuevo_potencial[i,j] = v_pico

            elif e == pra:  #No puede excitarse otra vez auneque el estimulo es fuerte(Si no, nos re morimos)
                nuevo_potencial[i,j] =  max(v_prr, V * (1 - t_pra_exp * dt)) 
                # --- PRA -> PRR ---
                if nuevo_potencial[i,j] <= v_prr: 
                    nuevo_estado[i,j] = prr

            elif e == prr:
                if I >= UP: #Si el estimulo es fuerte, puede volver a excitarse
                    nuevo_estado[i,j] = excitado
                    nuevo_potencial[i,j] = V #conserva el potencial actual al re-excitarse
                else:
                    nuevo_potencial[i,j] = max(v_reposo, V - t_rep * dt)
                    if nuevo_potencial[i,j] <= v_reposo: #cuando vuelve al potencial de reposo, pasa a reposo
                        nuevo_estado[i,j] = reposo
                        nuevo_potencial[i,j] = v_reposo   

    estado = nuevo_estado
    potencial = nuevo_potencial
    # Mostrar cada 20ms(para la grafica)
    if paso % 20 == 0:
        if paso>2000:
            break
        ax.clear()

        img = ax.imshow(
            estado,
            cmap=cmap,
            vmin=0,
            vmax=3
        )

        ax.set_title(f"Tiempo = {paso*dt:.0f} ms")

        plt.pause(0.01)
    potencial_celula.append(potencial[0, cols//2]) #registro del potencial de la celula estimulada a lo largo del tiempo
    ECG.append(potencial_total_vertical(potencial))

tiempos = np.arange(len(ECG)) * dt
#graficar_resultados(ECG, tiempos)
plt.figure(figsize=(10,4))
plt.plot(tiempos, potencial_celula, color='green')
plt.ylim(-100, 40)
plt.xlabel('Tiempo (ms)')
plt.ylabel('Potencial de membrana (mV)')
plt.title('Potencial de la célula marcapasos')
plt.grid(True)
plt.show()

