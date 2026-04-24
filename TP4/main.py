from modules.graficador import graficar_resultados, graficar_ingreso
from modules.metodos import metod_euler
import sympy as sp
import numpy as np

#Parametros generales del modelo --> la idea es que si el profe nos pide modificar algo solamente toquemos estos valores 
# y no cambiemos nada del codigo de la simulacion o de las graficas

parametros= {
    #valor de referencia de la glucosa en sangre
    'Gn': 100, #cantidad de glucosa por mililitro de sangre en ayunas

    #Variables endogenas del modelo
    'Is': sp.Symbol('Is'), #cantidad de insulina por mililitro de sangre 
    'Gs': sp.Symbol('Gs'), #cantidad de glucosa por mililitro de sangre 

    #Constantes del modelo
    'Koi': 0.005, #tasa de eliminacion de insulina por orina 
    'Kog': 0.0001, #tasa de eliminacion de glucosa por orina
    'kp': 0.1, #tasa de secrecion de insulina por parte del pancreas en respuesta a la glucosa
    'kh': 0.1, #tasa de secrecion de glucosa por parte del higado en respuesta a la glucosa (cuando esta por debajo de Gn)
    'kt': 0.0005, #tasa de eliminacion de glucosa por parte de los tejidos en respuesta a la insulina
    'qg': 0.1, #cantidad de glucosa ingerida por unidad de tiempo (puede ser 0 si no se ingiere nada)
    'qi': 0.1,  #cantidad de insulina inyectada por unidad de tiempo (puede ser 0 si no se inyecta nada)
    
    #Condicioens iniciales
    'Gs0': 100.0, #condicion inicial 
    'Is0': 5.0, #condicion inicial para insulina --> la guia dice em ayunas rango de 5 a 20 por eso cambie el 0
    
    #Simulacion
    'a': 0.0, # tiempo inicial
    'b': 1440, # tiempo de simulacion (1 dia en minutos)
    'dt': 1, #unidad de tiempo (minuto)

    'insulina_basal': 0.005, #cantidad de insulina basal para evitar que el modelo se vuelva inestable cuando Gs es muy bajo
    
     #ingesta de glucosa: (centro, amplitud, sigma) -> podemos cambiar horarios, intensidades y ancho de pulsos
     "ingestas": [
        (480, 80, 60),   # desayuno a las 8:00
        (780, 120, 60),  # almuerzo a las 13:00
        (1200, 100, 60),   # cena a las 20:00
     ],
     #inyeccion de insulina: (centro, amplitud, sigma) -> podemos cambiar horarios, intensidades y ancho de pulsos
     "inyecciones": [
        (470, 5, 20),   # bolo antes del desayuno
        (770, 7, 20),   # bolo antes del almuerzo
        (1190, 6, 20),   # bolo antes de la cena
     ],      
    }
    
# Parámetro para indicar si el paciente es diabético
diabetico = False  # Cambiar a True si el paciente es diabético
    
#-------------------------------------------------------
#Funciones auxiliares para la ingesta de glucosa y la inyeccion de insulina, modelizadas como "impulsos suaves" (gaussianas)
#----------------------------------------------------------

#Modelizamos la ingesta de glucosa como "impulsos suaves"
def pulsos (t,centro,amplitud,sigma):
    'Pulso sueva tipo gaussiano'
    'sirve para modelar tanto la ingesta de glucosa como la inyeccion de insulina, dependiendo de los parametros que le pasemos'
    'como un evento que sube y baja gradualmente, en vez de un salto brusco'
    return amplitud*np.exp(-(t-centro)**2/(2*sigma**2))   

def tren_de_pulsos(t, eventos):
    'Suma de pulsos para modelar multiples eventos a lo largo del dia'
    'eventos es una lista de tuplas (centro, amplitud, sigma) que definen cada pulso'
    return sum(pulsos(t, centro, amplitud, sigma) for centro, amplitud, sigma in eventos)  

def ingesta_glucosa(t, parametros=parametros):
    'Funcion de ingreso de glucosa '
    return tren_de_pulsos(t, parametros['ingestas'])

def inyeccion_insulina(t, parametros=parametros):
    'Funcion de ingreso de insulina'
    return tren_de_pulsos(t, parametros['inyecciones'])

#-----------------------------------------
#modelo compartimental
#------------------------------------------
def modelo_compartimental(t, y, parametros=parametros, diabetico=diabetico):
    Is, Gs = y
    Gn = parametros['Gn']
    Koi = parametros['Koi']
    Kog = parametros['Kog']
    kp = parametros['kp']
    kh = parametros['kh']
    kt = parametros['kt']
    qg = parametros['qg']
    qi = parametros['qi']

    #entradas externas
    Gin= ingesta_glucosa(t, parametros)

    #Agrego que si es diabetico se administre insulina externa segun la funcion de inyeccion_insulina,
    # sino solo la insulina basal para evitar que el modelo se vuelva inestable cuando Gs es muy bajo
    if diabetico:
        Iext = inyeccion_insulina(t, parametros) + parametros['insulina_basal'] #agregamos una pequeña cantidad de 
        #insulina basal para evitar que el modelo se vuelva inestable cuando Gs es muy bajo
    else:
        Iext = parametros['insulina_basal']  # no se administra insulina externa si no es diabetico

    #caso 1, glucosa por encima del nivel de referencia -> higado reduce su secrecion de glucosa
    if Gs > Gn: 
        dIs = kp*(Gs - Gn) - Koi*Is + qi*Iext
        dGs = qg*Gin - kt*Is*Gs - Kog*Gs
    
    #caso 2, glucosa por debajo del nivel de referencia -> higado aumenta su secrecion de glucosa
    else: 
        dIs = -Koi*Is + qi*Iext
        dGs = qg*Gin - kt*Is*Gs + kh*(Gn - Gs)
    return [dIs, dGs]

#----------------------------------------
#ejecutamos la simulacion
#------------------------------------------
simulacion = metod_euler(lambda t, y: modelo_compartimental(t, y, parametros, diabetico), 
                        parametros['a'], parametros['b'], [parametros['Is0'], 
                        parametros['Gs0']], parametros['dt'])
#graficamos resultados
graficar_resultados(simulacion, 'Modelo Compartimental Gs-Is')

#grafica de entradas externas
tiempos = [p[0] for p in simulacion]
Gin = [ingesta_glucosa(t, parametros) for t in tiempos]
if diabetico:
    Iext = [inyeccion_insulina(t, parametros) + parametros['insulina_basal'] for t in tiempos]
else:
    Iext = [parametros['insulina_basal'] for t in tiempos]
graficar_ingreso(tiempos, Gin, Iext)