import sympy as sp
import numpy as np


def metod_taylor_segundo_orden(f, dfdy, a, b, yi, h, dfdt=None):
    """Método de Taylor de orden 2 para y' = f(t, y).

    Parámetros:
    - f: función f(t, y) que devuelve y'.
    - dfdy: derivada de f con respecto a y (Jacobiano).
    - dfdt: opcional ∂f/∂t función de (t, y); si no se da, se asume vector de ceros.
    - a, b: intervalo de integración.
    - yi: valor inicial y(a).
    - h: paso numérico."""
    
    if dfdt is None:
        dfdt = lambda t, y: np.zeros_like(y)

    y = np.array(yi)
    t = a
    m = int((b - a) / h)
    lista = [(t, y.copy())]

    for i in range(m):
        f_val = np.array(f(t, y))
        jac = np.array(dfdy(t, y))
        derivada_segundo = np.array(dfdt(t, y)) + np.dot(jac, f_val)

        y = y + f_val * h + 0.5 * derivada_segundo * h**2
        y = np.maximum(y, 0)

        t = t + h
        lista.append((t, y.copy()))

    return lista
