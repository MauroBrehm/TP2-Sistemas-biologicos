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

def metod_ruger(f, a, b, yi, h, orden):
    m = int((b-a)/h)
    y = yi
    t = a
    if orden == 4:
        lista_RK4 = [(t, y)]
        for i in range(m):
            f1 = f(t, y)
            f2 = f(t+h/2, y+h/2*f1)
            f3 = f(t+h/2, y+h/2*f2)
            f4 = f(t+h, y+h*f3)
            y = y + (h*(f1+2*f2+2*f3+f4))/6
            if y < 0:
                y = 0
            t = t + h
            lista_RK4.append((t, y))
        return lista_RK4
    if orden == 2:
        ''' y = y+(a1*f1 + a2*f2)*h
        Donde, 
            f1 = f(t, y)
            f2 = f(t+p1*h, y+q11*f1*h)
        Aplicamos metodo del punto medio para dar valores a las constantes desconocidas
        a2 = 1 ; a1 = 0 ; p1 = q11 = 1/2 '''

        lista_RK2 = [(t,y)]
        for i in range(m):
            f_1 = f(t, y)
            f_2 = f(y+1/2*h, y+1/2*f_1*h)
            y = y + f_2*h
            if y < 0:
                y = 0
            t = t + h
            lista_RK2.append ((t, y))
        return lista_RK2
    
if __name__ == "__main__":
    pass