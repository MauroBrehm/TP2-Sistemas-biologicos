from modules.graficador import graficar_resultados
from modules.metodos import metod_euler
import sympy as sp
import numpy as np

Gn = 100 #cantidad de glucosa por mililitro de sangre en ayunas
"Variables endogenas del modelo" 
Is = sp.Symbol('Is') #cantidad de insulina por mililitro de sangre 
Gs = sp.Symbol('Gs') #cantidad de glucosa por mililitro de sangre 

#Parametros del modelo
Koi = 0.005 #tasa de eliminacion de insulina por orina 
Kog  = 0.0001#tasa de eliminacion de glucosa por orina
kp = 0.1 #tasa de secrecion de insulina por parte del pancreas en respuesta a la glucosa
kh = 0.1 #tasa de secrecion de glucosa por parte del higado en respuesta a la glucosa (cuando esta por debajo de Gn)
kt = 0.0005 #tasa de eliminacion de glucosa por parte de los tejidos en respuesta a la insulina
qg = 0.1 #cantidad de glucosa ingerida por unidad de tiempo (puede ser 0 si no se ingiere nada)
qi = 0.1 #cantidad de insulina inyectada por unidad de tiempo (puede ser 0 si no se inyecta nada)
Gs0 = 100 #condicion inicial 
b = 1440 # tiempo de simulacion (1 dia en minutos)
dt = 1 #unidad de tiempo (minuto)


#Modelizamos la ingesta de glucosa como "impulsos suaves"
def ingesta_glucosa(t):
    desayuno = 80*np.exp(-(t-480)**2/(2*60**2))   # 8:00
    almuerzo = 120*np.exp(-(t-780)**2/(2*60**2))  # 13:00
    cena = 100*np.exp(-(t-1200)**2/(2*60**2))     # 20:00
    return desayuno + almuerzo + cena 

#Modelizamos la inyeccion de insulina como "impulsos suaves"
def inyeccion_insulina(t):
    bolo_desayuno = 5*np.exp(-(t-470)**2/(2*20**2))
    bolo_almuerzo = 7*np.exp(-(t-770)**2/(2*20**2))
    bolo_cena = 6*np.exp(-(t-1190)**2/(2*20**2))
    return bolo_desayuno + bolo_almuerzo + bolo_cena
#simulacion
def modelización_compartimental_Gs_Is(t, y):
    Is, Gs = y
    Gin = ingesta_glucosa(t)
    Iext = inyeccion_insulina(t)+0.005 #agregamos una pequeña cantidad de insulina basal para evitar que el modelo se vuelva inestable cuando Gs es muy bajo

    if Gs > Gn:
        dIs = kp*(Gs - Gn) - Koi*Is + qi*Iext
        dGs = qg*Gin - kt*Is*Gs - Kog*Gs
    else:
        dIs = -Koi*Is + qi*Iext
        dGs = qg*Gin - kt*Is*Gs + kh*(Gn - Gs)
    return [dIs, dGs]
    
simulacion = metod_euler(modelización_compartimental_Gs_Is, 0, b, [Is, Gs], dt)