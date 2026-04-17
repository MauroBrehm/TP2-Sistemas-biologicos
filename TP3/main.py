
from modules.metodos import metod_euler
from modules.parametros import obtener_parametros
from modules.modelo import derivadas
from modules.graficador import graficar_resultados, graf_res

#Simulo sin la droga
datos=obtener_parametros()
#resuelvo el sistema de ecuaciones diferenciales con el metodo de euler
N = 1000
dt=0.01
a0=0
b0= N*dt

#valores iniciales
CD4_0 = datos["CD4N"]  # 1000
CD8_0 = datos["CD8N"]  # 550
V_0 = 0.001

def derivadas_para_retardo(t, y):
    droga_activa = t >= 2  # La droga se activa a partir del tiempo t=2 años
    return derivadas(y[0], y[1], y[2], datos, droga=droga_activa)

#Realizamos dos simulaciones:
#sin la droga, para ver la evolución natural del virus
simulacion_sin_droga= metod_euler(lambda t, y: derivadas(y[0], y[1], y[2], datos, droga=False), a0, b0, [CD4_0, CD8_0, V_0], dt)
#con la droga, para ver el efecto de la droga sobre la evolución del virus
simulacion_con_droga= metod_euler(derivadas_para_retardo, a0, b0, [CD4_0, CD8_0, V_0], dt)

# graficar_resultados(simulacion_sin_droga, 'Sin Droga')
# graficar_resultados(simulacion_con_droga, 'Con droga')

graf_res(simulacion_sin_droga, 'Sin Droga')
graf_res(simulacion_con_droga, 'Con Droga')
