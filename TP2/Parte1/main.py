import sympy as sp
import numpy as np
from modules.graficador import graficar_resultados, graficar_errores
from modules.metodos import metod_euler, metod_taylor_segundo_orden, metod_ruger
from modules.calculo_errores import calcular_error
# Parametros
t = sp.symbols('t')
P = sp.Function('P')
t_inicial = 0
t_final = 6 #años
P_inicial = 10
h = 1/12
h_2 = 2/12
h_3 = 1

# Definir la ecuación diferencial
eq = sp.Eq(sp.Derivative(P(t), t), (10)*P(t)*(1 - P(t)/1000))
#Resolver la EDO
sol = sp.dsolve(eq, P(t))

# Aplicar condición inicial P(0) = 10
cond_inic = {P(0): 10}
sol_exacta = sp.dsolve(eq,P(t), ics= cond_inic)
print(sol_exacta)

# Solución analítica evaluada
def P_exacta(t_val):
    return 1000 / (1 + 99 * sp.exp(-(10) * t_val))
t_vals = np.linspace(t_inicial, t_final, 1000)
real = [(t_val, float(P_exacta(t_val))) for t_val in t_vals]


#Ejercicio 2 
'''Metodo de EULER'''
def f(t, P):
        return (10)*P*(1-(P/1000))

res_euler_1 = metod_euler(f, t_inicial, t_final, P_inicial, h)
res_euler_2 = metod_euler (f, t_inicial, t_final, P_inicial, h_2)
res_euler_3 = metod_euler (f, t_inicial, t_final, P_inicial, h_3)
grafico = graficar_resultados([res_euler_1, res_euler_2, res_euler_3, real], ['Euler h=1 mes', 'Euler h=2 meses', 'Euler h=12 meses','Solución Exacta'], 'Euler')
error_E1 = calcular_error(res_euler_1, P_exacta)
error_E2 = calcular_error(res_euler_2, P_exacta)
error_E3 = calcular_error(res_euler_3, P_exacta)
graficar_errores([error_E1, error_E2, error_E3], ['h=1 mes', 'h=2 meses', 'h=12 meses'], 'Errores')

'''Metodo de TAYLOR de orden 2'''
res_taylor_1 = metod_taylor_segundo_orden(f, lambda t, P: 10*(1 - 2*P/1000), t_inicial, t_final, P_inicial, h)
res_taylor_2 = metod_taylor_segundo_orden(f, lambda t, P: 10*(1 - 2*P/1000), t_inicial, t_final, P_inicial, h_2)
res_taylor_3 = metod_taylor_segundo_orden(f, lambda t, P: 10*(1 - 2*P/1000), t_inicial, t_final, P_inicial, h_3)
grafico2 = graficar_resultados([res_taylor_1, res_taylor_2, res_taylor_3, real], ['Taylor h=1/12', 'Taylor h=2/12', 'Taylor h=1','Solución Exacta'], 'Taylor')
error_T1 = calcular_error(res_taylor_1, P_exacta)
error_T2 = calcular_error(res_taylor_2, P_exacta)
error_T3 = calcular_error(res_taylor_3, P_exacta)
graficar_errores([error_T1, error_T2, error_T3], ['h=1 mes', 'h=2 meses', 'h=12 meses'], 'Errores_Taylor')

'''Metodo de RUNGE_KUTTA orden 4'''
res_RK4_1 = metod_ruger(f, t_inicial, t_final, P_inicial, h, 4)
res_RK4_2 = metod_ruger(f,t_inicial, t_final, P_inicial, h_2, 4)
res_RK4_3 = metod_ruger(f, t_inicial, t_final, P_inicial, h_3, 4)
grafico_3 = graficar_resultados([res_RK4_1, res_RK4_2, res_RK4_3, real], ['RK4 h= 1/12', 'RK4 h=2/12', 'RK4 h=1', 'Solución exacta'], 'RK4')
error_RK41 = calcular_error(res_RK4_1, P_exacta)
error_RK42 = calcular_error(res_RK4_2, P_exacta)
error_RK43 = calcular_error(res_RK4_3, P_exacta)
graficar_errores([error_RK41, error_RK42, error_RK43], ['h=1 mes', 'h=2 meses', 'h=12 meses'], 'Errores_RK4')

'''Metodo de RUNGE-KUTTA orden 2'''
res_RK2_1 = metod_ruger(f, t_inicial, t_final, P_inicial, h, 2)
res_RK2_2 = metod_ruger(f, t_inicial, t_final, P_inicial, h_2, 2)
res_RK2_3 = metod_ruger(f, t_inicial, t_final, P_inicial, h_3, 2)
grafico_4 = graficar_resultados([res_RK2_1, res_RK2_2, res_RK2_3, real], ['RK2 h= 1/12', 'RK2 h=2/12', 'RK2 h=1', 'Solución exacta'], 'RK2')
error_RK21 = calcular_error(res_RK2_1, P_exacta)
error_RK22 = calcular_error(res_RK2_2, P_exacta)
error_RK23 = calcular_error(res_RK2_3, P_exacta)
graficar_errores([error_RK21, error_RK22, error_RK23], ['h=1 mes', 'h=2 meses', 'h=12 meses'], 'Errores RK2')

'''Grafico de paso h=1 en distintos metodos''' #era solo comparar con los distintos metodos el paso h=1 (año)
graficar_resultados([res_euler_3, res_RK2_3, res_RK4_3, res_taylor_3], ['Euler', 'RK2', 'RK4', 'Taylor'], 'paso1')