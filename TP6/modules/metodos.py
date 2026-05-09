import numpy as np

# Método de Euler para resolver sistemas de EDOs
def metod_euler(f , a, b, yi, h):
    y = np.array(yi) #valor anterior
    t = a  # tiempo inicial
    m = int((b-a)/h)  # número de pasos
    lista = [(t, y.copy())]  #se usa copy para guardar una copia del valor en cada paso, 
                             #sino hacemos esto los valores se sobreescribirían y todas las graficas serian lineas
    for i in range(m):
        y = y + np.array(f(t, y)) * h 
        t = t + h 
        lista.append((t, y.copy()))  
    return lista

