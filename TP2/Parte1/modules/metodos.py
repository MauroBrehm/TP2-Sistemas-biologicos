import sympy as sp
import numpy as np

def metod_euler(f , a, b, yi, h):
    y = yi #valor anterior
    t = a  # tiempo inicial
    m = int((b-a)/h)  # número de pasos
    lista = [(t, y)]  
    for i in range(m):
        y = y + f(t, y) * h 
        if y < 0:
            y = 0
        t = t + h 
        lista.append((t, y))  
    return lista

def metod_taylor_segundo_orden(f, dfdy, a, b, yi, h, dfdt=None):
    """Método de Taylor de orden 2 para y' = f(t, y).

    

    Parámetros:
    - f: función f(t, y) que devuelve y'.
    - dfdy: derivada de f con respecto a y.
    - dfdt: opcional ∂f/∂t función de (t, y); si no se da, se asume 0.
    - a, b: intervalo de integración.
    - yi: valor inicial y(a).
    - h: paso numérico."""
    
    if dfdt is None:
        dfdt = lambda t, y: 0

    y = yi
    t = a
    m = int((b - a) / h)
    lista = [(t, y)]

    for i in range(m):
        f_val = f(t, y)
        derivada_segundo = dfdt(t, y) + dfdy(t, y) * f_val

        y = y + f_val * h + 0.5 * derivada_segundo * h**2
        if y < 0:
            y = 0

        t = t + h
        lista.append((t, y))

    return lista




def metod_ruger():
    pass


if __name__ == "__main__":
    pass