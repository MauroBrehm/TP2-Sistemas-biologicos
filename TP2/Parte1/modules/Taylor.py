import sympy as sp
from modules.funcion import P




def serieTaylor(f, a, n):
    x = sp.symbols('x', real=True)
    serie = 0
    for k in range(n+1):
        derivada = sp.diff(f, x, k).subs(x, a)
        serie += derivada / sp.factorial(k) * (x - a)**k
    return serie

