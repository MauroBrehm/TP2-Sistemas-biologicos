from modules.metodo import metod_taylor_segundo_orden

import matplotlib.pyplot as plt
#Definimos la funcion que es el Modelo SIR
def f(t,N):
    S, I, R = N
    beta = 0.3
    gamma = 0.1
    N_total = 1000
    dSdt = -beta * S * I / N_total
    dIdt = beta * S * I / N_total - gamma * I
    dRdt = gamma * I
    return [dSdt, dIdt, dRdt]

#Definimos el Jacobiano para el método de Taylor
def jacobian(t, N):
    S, I, R = N
    beta = 0.3
    gamma = 0.1
    N_total = 1000
    return [
        [-beta * I / N_total, -beta * S / N_total, 0],
        [beta * I / N_total, beta * S / N_total - gamma, 0],
        [0, gamma, 0]
    ]

#Implementamos Taylor de orden 2 para el modelo SIR 
#paso=0.5 y N=1000
#Condiciones iniciales: S(0)=990, I(0)=10, R(0)=0 durante 60 dias
a = 0
b = 60
h = 0.5
N = 1000
S0 = 990
I0 = 10
R0 = 0
resultados = metod_taylor_segundo_orden(f, jacobian, a, b, [S0, I0, R0], h)
#Imprimimos los resultados
for t, N in resultados:
    print(f"t={t:.1f}, S={N[0]:.2f}, I={N[1]:.2f}, R={N[2]:.2f}")
#Verificar que el numero basico de reproduccion es R0=beta/gamma=3.0
R0 = 0.3 / 0.1
print(f"El número básico de reproducción R0 es: {R0:.2f}")
#Mostramos I(t) en 30 puntos equispaciados entre 0 y 60


t_values = [t for t, _ in resultados]
I_values = [N[1] for _, N in resultados]

plt.plot(t_values, I_values, label='I(t)')
plt.xlabel('t')
plt.ylabel('I')
plt.title('Evolución de la población infectada')
plt.legend()
plt.show()

#Implementación en 16 bits de la ecuacion:𝑝 = 𝑝min + (𝑝𝑚𝑎𝑥 −𝑝𝑚𝑖𝑛)∗ 𝑑𝑒𝑐𝑖𝑚𝑎𝑙(𝑏₁𝑏₂ … 𝑏. . 𝐿)/(2**L- 1) 
#Verificar!! que el cromosoma de todos 0 da Beta_min y Gama_min y el cromosoma de todos 1 da el Beta_max y Gama_max

def decode_chromosome(chromosome, p_min, p_max):
    L = len(chromosome)
    decimal_value = sum(bit * (2 ** (L - 1 - i)) for i, bit in enumerate(chromosome))
    p = p_min + (p_max - p_min) * decimal_value / (2**L - 1)
    return p
