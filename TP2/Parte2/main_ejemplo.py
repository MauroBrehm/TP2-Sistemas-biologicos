from modules.metodo import metod_taylor_segundo_orden
import random
import matplotlib.pyplot as plt
from modules.operadores_geneticos import encode_beta_gamma, decode_beta_gamma, calcular_ecm

#Definimos la funcion que es el Modelo SIR
beta = 0.3
gamma = 0.1
def f(t, N):
    '''Modelo SIR'''
    S, I, R = N
    # beta = 0.3
    # gamma = 0.1
    N_total = 1000
    dSdt = -beta * S * I / N_total
    dIdt = beta * S * I / N_total - gamma * I
    dRdt = gamma * I
    return [dSdt, dIdt, dRdt]

#Definimos el Jacobiano para el método de Taylor
def jacobian(t, N):
    S, I, R = N
    # beta = 0.3
    # gamma = 0.1
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

print ('-'*50)
print("Implementacion del AG - Codificacion y Decodificación ")
#Implementación en 16 bits de la ecuacion:𝑝 = 𝑝min + (𝑝𝑚𝑎𝑥 −𝑝𝑚𝑖𝑛)∗ 𝑑𝑒𝑐𝑖𝑚𝑎𝑙(𝑏₁𝑏₂ … 𝑏. . 𝐿)/(2**L- 1) 

# rangos de busqueda de beta y gamma
beta_min = 0.05
beta_max = 1
gamma_min = 0.01
gamma_max = 0.5

# Verificación decodificación
cromosoma_ceros= [0] * 16
beta_dec, gamma_dec = decode_beta_gamma(cromosoma_ceros, beta_min, beta_max, gamma_min, gamma_max)
print(f"Cromosoma todos ceros: beta={beta_dec}, gamma={gamma_dec}")

cromosoma_unos = [1] * 16
beta_dec, gamma_dec = decode_beta_gamma(cromosoma_unos, beta_min, beta_max, gamma_min, gamma_max)
print(f"Cromosoma todos unos: beta={beta_dec}, gamma={gamma_dec}")

# Generar uno random
cromosoma_random = [random.randint(0, 1) for _ in range(16)]
beta_dec, gamma_dec = decode_beta_gamma(cromosoma_random, beta_min, beta_max, gamma_min, gamma_max)
print(f"Cromosoma random: {cromosoma_random}")
print(f"Decodificado cromosoma random: beta={beta_dec}, gamma={gamma_dec}")

#Implementacion de funcion fitness
print('-'*50)
print ('A continuacion se muestra la implementacion del fitness')

# Datos observados: usamos la simulación con beta=0.3, gamma=0.1 como referencia
I_obs = [N[1] for _, N in resultados]


print('Se integra el modelo SIR para cada individuo (β, γ)')
individuos = {
    'ceros': decode_beta_gamma(cromosoma_ceros, beta_min, beta_max, gamma_min, gamma_max),
    'unos': decode_beta_gamma(cromosoma_unos, beta_min, beta_max, gamma_min, gamma_max),
    'verdadero': (0.3, 0.1),  # Valores verdaderos para verificar fitness=1
    'random': decode_beta_gamma(cromosoma_random, beta_min, beta_max, gamma_min, gamma_max),
}

for nombre, (beta_val, gamma_val) in individuos.items():
    beta, gamma = beta_val, gamma_val
    res_individuo = metod_taylor_segundo_orden(f, jacobian, a, b, [S0, I0, R0], h)
    I_simulado = [N[1] for _, N in res_individuo] # datos de curva infectadors integrando SIR
    ecm = calcular_ecm(I_simulado, I_obs)
    fitness = 1/(1+ecm)
    
    print(f"Individuo {nombre}: beta={beta_val:.2f}, gamma={gamma_val:.2f}, ECM={ecm:.2f}, fitness={fitness}")
    
