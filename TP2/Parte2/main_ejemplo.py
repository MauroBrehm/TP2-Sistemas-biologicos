from modules.metodo import metod_taylor_segundo_orden
import random
import matplotlib.pyplot as plt
from modules.operadores_geneticos import encode_beta_gamma, decode_beta_gamma, calcular_ecm
from modules.operadores_geneticos import seleccion_ruleta, mutacion, cruzamiento

#Esta funcion va a generar el modelo SIR para valores especificos de beta y gamma
#se usa para que cada individuo del AG tenga su propio modelo
def crear_modelo_sir(beta, gamma):
    def f(t, N): #sistema EDO
        S, I, R = N
        N_total = 1000
        dSdt = -beta * S * I / N_total
        dIdt = beta * S * I / N_total - gamma * I
        dRdt = gamma * I
        return [dSdt, dIdt, dRdt]

#Definimos el Jacobiano para el método de Taylor
    def jacobian(t, N):
        S, I, R = N  
        N_total = 1000
        return [
            [-beta * I / N_total, -beta * S / N_total, 0],
            [beta * I / N_total, beta * S / N_total - gamma, 0],
            [0, gamma, 0]
        ]
    return f, jacobian

#Implementamos Taylor de orden 2 para el modelo SIR 
#intervalo de simulacion(dias)
a = 0
b = 60
#paso
h = 0.5
#poblacion total
N = 1000
#condiciones iniciales
S0 = 990
I0 = 10
R0 = 0
#generamos datos "reales" usando los valores verdaderos de beta y gamma, van a ser los que el AG intente reproducir
f_real,jac_real = crear_modelo_sir(0.3, 0.1)

resultados = metod_taylor_segundo_orden(f_real, jac_real, a, b, [S0, I0, R0], h)

#Imprimimos los resultados
for t, estado in resultados:
    print(f"t={t:.1f}, S={estado[0]:.2f}, I={estado[1]:.2f}, R={estado[2]:.2f}")

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

# Datos observados de infectados: usamos la simulación con beta=0.3, gamma=0.1 como referencia
#pedia solo 30 pero nosotros tenemos como 120 puntos por eso lo cambie pero igual lo dejo comentado por si las dudas
#I_obs = [N[1] for _, N in resultados]
indices=[int(i*(len(resultados)-1)/29) for i in range(30)]
I_obs=[resultados[i][1][1] for i in indices]

#verificamos fitness
beta_test, gamma_test = 0.3, 0.1
f_test, jac_test = crear_modelo_sir(beta_test, gamma_test)

resultados_test = metod_taylor_segundo_orden(f_test, jac_test, a, b, [S0, I0, R0], h)
I_test_total = [estado[1] for _, estado in resultados_test]
I_test = [I_test_total[i] for i in indices]

ecm_test = calcular_ecm(I_test, I_obs)
fitness_test = 1 / (1 + ecm_test)
print(f"Fitness con valores verdaderos (beta=0.3, gamma=0.1): {fitness_test:.4f}")

print('Se integra el modelo SIR para cada individuo (β, γ)')

#parametros del AG
#estos valores los saque del final del TP no se si esta bien
n_poblacion = 100 #cant de individuos
L=16 #longitud del cromosoma (bits para beta + bits para gamma)
generaciones = 50 #cant de generaciones del AG (iteraciones)

#poblacion inicial (aleatoria)
poblacion=[
    [random.randint(0, 1) for _ in range(L)]
    for _ in range(n_poblacion)
]

#funcion de fitness
#evalua que tan bueno es cada individuo de la poblacion
def evaluar_poblacion(poblacion):
    fitnesses = []
    for cromosoma in poblacion:
        #decodificamos el cromosoma para obtener beta y gamma
        beta_val, gamma_val = decode_beta_gamma(cromosoma, beta_min, beta_max, gamma_min, gamma_max)
        #creamos el modelo SIR para esos valores
        f_individuo, jac_indiviuo = crear_modelo_sir(beta_val, gamma_val)
        #simulamos el modelo con Taylor de orden 2
        resultados_sim = metod_taylor_segundo_orden(f_individuo, jac_indiviuo, a, b, [S0, I0, R0], h)
        #sacamos curva de infectados simulados
        I_sim_totsl = [estado[1] for _, estado in resultados_sim]
        I_sim = [I_sim_totsl[i] for i in indices] #solo los 30 puntos que tenemos de observados
        #calculamos el ECM entre la curva simulada y la observada
        ecm = calcular_ecm(I_sim, I_obs)
        fitness = 1 / (1 + ecm)  #Evitamos división por cero
        fitnesses.append(fitness)
    return fitnesses

#listas para guardar el mejor fitness y el promedio de cada generación
mejores=[]
promedios=[]

#bucle evolutivo del AG
for gen in range(generaciones):
    print(f"Generación {gen + 1}")
    fitnesses = evaluar_poblacion(poblacion)
    mejores.append(max(fitnesses))
    promedios.append(sum(fitnesses) / len(fitnesses))

    #elitismo
    mejor_indice = fitnesses.index(max(fitnesses))
    mejor_individuo = poblacion[mejor_indice][:]
    
    #seccionamos usando ruleta, aplicamos cruzamiento y mutacion para generar la nueva poblacion
    seleccionados = seleccion_ruleta(poblacion, fitnesses, n_poblacion - 1)
    nueva_poblacion=[]
    for i in range (0, n_poblacion, 2):
        padre1=seleccionados[i]
        padre2=seleccionados[i+1]
        #cruzamiento de un punto
        hijo1, hijo2 = cruzamiento(padre1, padre2)
        #mutacion de los hijos
        nueva_poblacion.append(mutacion(hijo1))
        nueva_poblacion.append(mutacion(hijo2))
    #aseguramos que el mejor individuo de la generación anterior se mantenga en la nueva población (elitismo)
    nueva_poblacion[0]=mejor_individuo
    #actualizamos la poblacion para la siguiente generacion
    poblacion=nueva_poblacion

#resultados finales
fitnesses_final = evaluar_poblacion(poblacion)

#obtener el mejor individuo de final
mejor_indice_final = fitnesses_final.index(max(fitnesses_final))
mejor=poblacion[mejor_indice_final]

#decodificamos el mejor individuo para obtener los valores de beta y gamma
beta_final, gamma_final = decode_beta_gamma(mejor, beta_min, beta_max, gamma_min, gamma_max)
print("\nMejor solución encontrada:")
print(f"beta: {beta_final}, gamma: {gamma_final}")
print(f"R0 = {beta_final / gamma_final:.2f}")

#error relativo porcentual
error_beta = abs(beta_final - 0.3) / 0.3 * 100
error_gamma = abs(gamma_final - 0.1) / 0.1 * 100
print(f"Error relativo porcentual en beta: {error_beta:.2f}%")
print(f"Error relativo porcentual en gamma: {error_gamma:.2f}%")

#grafico del AG
plt.plot(mejores, label='Mejor fitness')
plt.plot(promedios, label='Fitness promedio')
plt.xlabel('Generación')
plt.ylabel('Fitness')
plt.legend()
plt.title("Evolución del algoritmo genético")
plt.show()

#simulacion final con beta* y gamma*
f_final, jac_final = crear_modelo_sir(beta_final, gamma_final)
resultados_final = metod_taylor_segundo_orden(f_final, jac_final, a, b, [S0, I0, R0], h)
t_final = [t for t, _ in resultados_final]
S_final = [estado[0] for _, estado in resultados_final]
I_final = [estado[1] for _, estado in resultados_final]
R_final = [estado[2] for _, estado in resultados_final]

#grafico comparativo
plt.plot(t_final, S_final, label='S(t)')
plt.plot(t_final, I_final, label='I(t) estimado')
plt.plot(t_final, R_final, label='R(t)')

#graficamos I(t) observado para comparar
t_obs = [t_final[i] for i in indices]
plt.scatter(t_obs, I_obs, color='pink', label='I(t) observado')

plt.xlabel('tiempo')
plt.ylabel('población')
plt.legend()
plt.title('Comparación entre modelo SIR estimado vs datos observados')
plt.show()