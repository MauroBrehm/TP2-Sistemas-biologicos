import sympy as sp
import numpy as np
from modules.graficador import graficar_resultados, graficar_errores, graficar_comparacion
from modules.metodos import metod_euler, metod_taylor_segundo_orden, metod_ruger, metod_euler32,metod_taylor_segundo_orden32,metod_ruger32
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

#parametros para float 32
t_inicial32 = np.float32(0)
t_final32 = np.float32(6) #años
P_inicial32 = np.float32(10)
h32 = np.float32(1/12)
h_2_32 = np.float32(2/12)
h_3_32 = np.float32(1)

# Definir la ecuación diferencial
eq = sp.Eq(sp.Derivative(P(t), t), (10)*P(t)*(1 - P(t)/1000))
#Resolver la EDO
sol = sp.dsolve(eq, P(t))

#definir la ecuacion en float 32
eq32= sp.Eq(sp.Derivative(P(t), t), (np.float32(10))*P(t)*(np.float32(1) - P(t)/np.float32(1000)))
sol32 = sp.dsolve(eq32, P(t))


# Aplicar condición inicial P(0) = 10
cond_inic = {P(0): 10}
sol_exacta = sp.dsolve(eq,P(t), ics= cond_inic)
print(sol_exacta)

#Aplicar condiciones iniciales a la solución en float 32
cond_inic32 = {P(0): np.float32(10)}
sol_exacta32 = sp.dsolve(eq32, P(t), ics=cond_inic32)


# Solución analítica evaluada
def P_exacta(t_val):
    return 1000 / (1 + 99 * sp.exp(-(10) * t_val))
t_vals = np.linspace(t_inicial, t_final, 1000)
real = [(t_val, float(P_exacta(t_val))) for t_val in t_vals]

#Solución analitica evaluada en float 32
def P_exacta32(t_val):
    return np.float32(1000) / (np.float32(1) + np.float32(99) * np.exp(-np.float32(10) * t_val))
t_vals32 = np.linspace(t_inicial32, t_final32, 1000, dtype=np.float32)
real32 = [(t_val, np.float32(P_exacta32(t_val))) for t_val in t_vals32]

#Ejercicio 2 
'''Metodo de EULER'''
def f(t, P):
        return (10)*P*(1-(P/1000))

def f32(t, P):
        return (np.float32(10))*P*(np.float32(1)-(P/np.float32(1000)))

res_euler_1 = metod_euler(f, t_inicial, t_final, P_inicial, h)
res_euler_2 = metod_euler (f, t_inicial, t_final, P_inicial, h_2)
res_euler_3 = metod_euler (f, t_inicial, t_final, P_inicial, h_3)
grafico = graficar_resultados([res_euler_1, res_euler_2, res_euler_3, real], ['Euler h=1 mes', 'Euler h=2 meses', 'Euler h=12 meses','Solución Exacta'], 'Euler')
error_E1 = calcular_error(res_euler_1, P_exacta)
error_E2 = calcular_error(res_euler_2, P_exacta)
error_E3 = calcular_error(res_euler_3, P_exacta)
graficar_errores([error_E1, error_E2, error_E3], ['h=1 mes', 'h=2 meses', 'h=12 meses'], 'Errores')

#Ahora con float 32 y luego comparamos los errores
res_euler_1_32 = metod_euler32(f32, t_inicial32, t_final32, P_inicial32, h32)
res_euler_2_32 = metod_euler32(f32, t_inicial32, t_final32, P_inicial32, h_2_32)
res_euler_3_32 = metod_euler32(f32, t_inicial32, t_final32, P_inicial32, h_3_32)
#print(type(res_euler_1_32[0][1])) #verificar que es float 32
error_E1_32 = calcular_error(res_euler_1_32, P_exacta32)
error_E2_32 = calcular_error(res_euler_2_32, P_exacta32)
error_E3_32 = calcular_error(res_euler_3_32, P_exacta32)
#print(type(error_E1_32[0][1])) #verificar que es float 32
diferencia_error_Euler_h1 = [(t, abs(e32 - e64)) for ((t, e32), (_, e64)) in zip(error_E1_32, error_E1)]
diferencia_error_Euler_h2 = [(t, abs(e32 - e64)) for ((t, e32), (_, e64)) in zip(error_E2_32, error_E2)]
diferencia_error_Euler_h3 = [(t, abs(e32 - e64)) for ((t, e32), (_, e64)) in zip(error_E3_32, error_E3)]
graficar_comparacion([diferencia_error_Euler_h1, diferencia_error_Euler_h2, diferencia_error_Euler_h3], ['h=1 mes', 'h=2 meses', 'h=12 meses'], 'Diferencia errores Euler')
'''Metodo de TAYLOR de orden 2'''
res_taylor_1 = metod_taylor_segundo_orden(f, lambda t, P: 10*(1 - 2*P/1000), t_inicial, t_final, P_inicial, h)
res_taylor_2 = metod_taylor_segundo_orden(f, lambda t, P: 10*(1 - 2*P/1000), t_inicial, t_final, P_inicial, h_2)
res_taylor_3 = metod_taylor_segundo_orden(f, lambda t, P: 10*(1 - 2*P/1000), t_inicial, t_final, P_inicial, h_3)
grafico2 = graficar_resultados([res_taylor_1, res_taylor_2, res_taylor_3, real], ['Taylor h=1/12', 'Taylor h=2/12', 'Taylor h=1','Solución Exacta'], 'Taylor')
error_T1 = calcular_error(res_taylor_1, P_exacta)
error_T2 = calcular_error(res_taylor_2, P_exacta)
error_T3 = calcular_error(res_taylor_3, P_exacta)
graficar_errores([error_T1, error_T2, error_T3], ['h=1 mes', 'h=2 meses', 'h=12 meses'], 'Errores_Taylor')

"Metodo de Taylor de orden 2 con float 32"
res_taylor_1_32 = metod_taylor_segundo_orden32(f32, lambda t, P: 10*(1 - 2*P/1000), t_inicial32, t_final32, P_inicial32, h32)
res_taylor_2_32 = metod_taylor_segundo_orden32(f32, lambda t, P: 10*(1 - 2*P/1000), t_inicial32, t_final32, P_inicial32, h_2_32)
res_taylor_3_32 = metod_taylor_segundo_orden32(f32, lambda t, P: 10*(1 - 2*P/1000), t_inicial32, t_final32, P_inicial32, h_3_32)
#print
eror_T1_32 = calcular_error(res_taylor_1_32, P_exacta32)
eror_T2_32 = calcular_error(res_taylor_2_32, P_exacta32)
eror_T3_32 = calcular_error(res_taylor_3_32, P_exacta32)
#print(type(eror_T1_32[0][1])) #verificar que es float 32
diferencia_error_Taylor_h1 = [(t, abs(e32 - e64)) for ((t, e32), (_, e64)) in zip(eror_T1_32, error_T1)]
diferencia_error_Taylor_h2 = [(t, abs(e32 - e64)) for ((t, e32), (_, e64)) in zip(eror_T2_32, error_T2)]
diferencia_error_Taylor_h3 = [(t, abs(e32 - e64)) for ((t, e32), (_, e64)) in zip(eror_T3_32, error_T3)]
graficar_comparacion([diferencia_error_Taylor_h1, diferencia_error_Taylor_h2, diferencia_error_Taylor_h3], ['h=1 mes', 'h=2 meses', 'h=12 meses'], 'Diferencia errores Taylor')

'''Metodo de RUNGE_KUTTA orden 4'''
res_RK4_1 = metod_ruger(f, t_inicial, t_final, P_inicial, h, 4)
res_RK4_2 = metod_ruger(f,t_inicial, t_final, P_inicial, h_2, 4)
res_RK4_3 = metod_ruger(f, t_inicial, t_final, P_inicial, h_3, 4)
grafico_3 = graficar_resultados([res_RK4_1, res_RK4_2, res_RK4_3, real], ['RK4 h= 1/12', 'RK4 h=2/12', 'RK4 h=1', 'Solución exacta'], 'RK4')
error_RK41 = calcular_error(res_RK4_1, P_exacta)
error_RK42 = calcular_error(res_RK4_2, P_exacta)
error_RK43 = calcular_error(res_RK4_3, P_exacta)
graficar_errores([error_RK41, error_RK42, error_RK43], ['h=1 mes', 'h=2 meses', 'h=12 meses'], 'Errores_RK4')

"Metodo de RUNGE-KUTTA orden 4 con float 32"
res_RK4_1_32 = metod_ruger32(f32, t_inicial32, t_final32, P_inicial32, h32, 4)
res_RK4_2_32 = metod_ruger32(f32, t_inicial32, t_final32, P_inicial32, h_2_32, 4)
res_RK4_3_32 = metod_ruger32(f32, t_inicial32, t_final32, P_inicial32, h_3_32, 4)
#print(type(res_RK4_1_32[0][1])) #verificar que es float 32
error_RK41_32 = calcular_error(res_RK4_1_32, P_exacta32)
error_RK42_32 = calcular_error(res_RK4_2_32, P_exacta32)
error_RK43_32 = calcular_error(res_RK4_3_32, P_exacta32)
#print(type(error_RK41_32[0][1])) #verificar que es float 32
diferencia_error_RK4_h1 = [(t, abs(e32 - e64)) for ((t, e32), (_, e64)) in zip(error_RK41_32, error_RK41)]
diferencia_error_RK4_h2 = [(t, abs(e32 - e64)) for ((t, e32), (_, e64)) in zip(error_RK42_32, error_RK42)]
diferencia_error_RK4_h3 = [(t, abs(e32 - e64)) for ((t, e32), (_, e64)) in zip(error_RK43_32, error_RK43)]
graficar_comparacion([diferencia_error_RK4_h1, diferencia_error_RK4_h2, diferencia_error_RK4_h3], ['h=1 mes', 'h=2 meses', 'h=12 meses'], 'Diferencia errores RK4')

'''Metodo de RUNGE-KUTTA orden 2'''
res_RK2_1 = metod_ruger(f, t_inicial, t_final, P_inicial, h, 2)
res_RK2_2 = metod_ruger(f, t_inicial, t_final, P_inicial, h_2, 2)
res_RK2_3 = metod_ruger(f, t_inicial, t_final, P_inicial, h_3, 2)
grafico_4 = graficar_resultados([res_RK2_1, res_RK2_2, res_RK2_3, real], ['RK2 h= 1/12', 'RK2 h=2/12', 'RK2 h=1', 'Solución exacta'], 'RK2')
error_RK21 = calcular_error(res_RK2_1, P_exacta)
error_RK22 = calcular_error(res_RK2_2, P_exacta)
error_RK23 = calcular_error(res_RK2_3, P_exacta)
graficar_errores([error_RK21, error_RK22, error_RK23], ['h=1 mes', 'h=2 meses', 'h=12 meses'], 'Errores RK2')

"Metodo de RUNGE-KUTTA orden 2 con float 32"
res_RK2_1_32 = metod_ruger32(f32, t_inicial32, t_final32, P_inicial32, h32, 2)
res_RK2_2_32 = metod_ruger32(f32, t_inicial32, t_final32, P_inicial32, h_2_32, 2)
res_RK2_3_32 = metod_ruger32(f32, t_inicial32, t_final32, P_inicial32, h_3_32, 2)
#print(type(res_RK2_1_32[0][1])) #verificar que es float 32
error_RK21_32 = calcular_error(res_RK2_1_32, P_exacta32)
error_RK22_32 = calcular_error(res_RK2_2_32, P_exacta32)
error_RK23_32 = calcular_error(res_RK2_3_32, P_exacta32)
#print(type(error_RK21_32[0][1])) #verificar que es float 32
diferencia_error_RK2_h1 = [(t, abs(e32 - e64)) for ((t, e32), (_, e64)) in zip(error_RK21_32, error_RK21)]
diferencia_error_RK2_h2 = [(t, abs(e32 - e64)) for ((t, e32), (_, e64)) in zip(error_RK22_32, error_RK22)]
diferencia_error_RK2_h3 = [(t, abs(e32 - e64)) for ((t, e32), (_, e64)) in zip(error_RK23_32, error_RK23)]
graficar_comparacion([diferencia_error_RK2_h1, diferencia_error_RK2_h2, diferencia_error_RK2_h3], ['h=1 mes', 'h=2 meses', 'h=12 meses'], 'Diferencia errores RK2')

'''Grafico de paso h=1 en distintos metodos''' #era solo comparar con los distintos metodos el paso h=1 (año)
graficar_resultados([res_euler_3, res_RK2_3, res_RK4_3, res_taylor_3], ['Euler', 'RK2', 'RK4', 'Taylor'], 'paso1')