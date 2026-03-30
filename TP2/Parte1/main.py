import sympy as sp
import numpy as np
from modules.graficador import grafico_euler,grafico_taylor
from modules.metodos import metod_euler, metod_taylor_segundo_orden
from modules.calculo_errores import calcular_error
# Parametros
t = sp.symbols('t')
P = sp.Function('P')
t_inicial = 0
t_final = 60
P_inicial = 10
h = 1/12
h_2 = 2/12
h_3 = 1

# Definir la ecuación diferencial
eq = sp.Eq(sp.Derivative(P(t), t), 10*P(t)*(1 - P(t)/1000))
# Resolver la EDO
sol = sp.dsolve(eq, P(t))
# Aplicar condición inicial P(0) = 10
cond_inic = {P(0): 10}
sol_exacta = sp.dsolve(eq,P(t), ics= cond_inic)
print(sol_exacta)

# Solución analítica evaluada
def P_exacta(t_val):
    return 1000 / (1 + 99 * sp.exp(-10 * t_val))
t_vals = np.linspace(t_inicial, t_final, 1000)
real = [(t_val, float(P_exacta(t_val))) for t_val in t_vals]

#Ejercicio 2 
'''Metodo de EULER'''
def f(t, P):
        return 10*P*(1-(P/1000))

res_euler_1 = metod_euler(f, t_inicial, t_final, P_inicial, h)
res_euler_2 = metod_euler (f, t_inicial, t_final, P_inicial, h_2)
res_euler_3 = metod_euler (f, t_inicial, t_final, P_inicial, h_3)
grafico = grafico_euler([res_euler_1, res_euler_2, res_euler_3, real], ['Euler h=1/12', 'Euler h=2/12', 'Euler h=1','Solución Exacta'])

"Metodo de taylor de orden 2"
res_taylor_1 = metod_taylor_segundo_orden(f, lambda t, P: 10*(1 - 2*P/1000), t_inicial, t_final, P_inicial, h)
res_taylor_2 = metod_taylor_segundo_orden(f, lambda t, P: 10*(1 - 2*P/1000), t_inicial, t_final, P_inicial, h_2)
res_taylor_3 = metod_taylor_segundo_orden(f, lambda t, P: 10*(1 - 2*P/1000), t_inicial, t_final, P_inicial, h_3)
grafico2 = grafico_taylor([res_taylor_1, res_taylor_2, res_taylor_3, real], ['Taylor h=1/12', 'Taylor h=2/12', 'Taylor h=1','Solución Exacta'])


