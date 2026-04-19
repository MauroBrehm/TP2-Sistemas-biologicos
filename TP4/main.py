from modules.graficador import graficar_resultados
from modules.metodos import metod_euler


Gn = 100 #cantidad de glucosa por mililitro de sangre en ayunas 
Gs = [Gn]  #no siento q este bien pero sino no sabia como hacer
Is = [5] #cantidad de insulina por mililitro de sangre en ayunas (no siento q este bien)

#Parametros del modelo
Kio = 0.005
Kog  = 0.0001
kp = 0.1
kh = 0.1
kt = 0.0005
qg = 0.1
qi = 0.1
Gs0 = 100 #condicion inicial 
b = 1440 # tiempo de simulacion (1 dia en minutos)
dt = 1 #unidad de tiempo (minuto)
'''No siento que esto este bien asi q se fijan jeje''' #lo puse asi para que vean 
Iin = 0.1 #cantidad de insulina inyectada por unidad de tiempo
Gin = 0.1 #cantidad de glucosa ingerida por unidad de tiempo


#simulacion
def derivadas (t, y):
    Is, Gs = y
    if Gs > Gn:
        dIs = kp*(Gs - Gn) - Kio*Is + qi*Iin
        dGs = qg*Gin - kt*Is*Gs - Kog*Gs
    else:
        dIs = -Kio*Is + qi*Iin
        dGs = qg*Gin - kt*Is*Gs + kh*(Gn - Gs)
    return [dIs, dGs]
    
simulacion = metod_euler(derivadas, 0, b, [Is, Gs], dt)