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

def serieTaylor(f, a, n):
    x = sp.symbols('x', real=True)
    serie = 0
    for k in range(n+1):
        derivada = sp.diff(f, x, k).subs(x, a)
        serie += derivada / sp.factorial(k) * (x - a)**k
    return serie


def metod_ruger():
    pass