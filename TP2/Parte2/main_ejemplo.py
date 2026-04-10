
import random
from modules.metodo import metod_taylor_segundo_orden, crear_modelo_sir
from modules.operadores_geneticos import decode_beta_gamma, calcular_ecm, evaluar_poblacion, seleccion_ruleta, mutacion, cruzamiento, ejecutar_ag_con_rangos
from modules.graficador import graficar_modelo_AG, grafico_comparacion_curvas, graficar_vector_I

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

graficar_vector_I(t_values, I_values)

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
#listas para guardar el mejor fitness y el promedio de cada generación
mejores=[]
promedios=[]
print('Se crean las generaciones')
#bucle evolutivo del AG
for gen in range(generaciones):
    # print(f"Generación {gen + 1}")
    fitnesses = evaluar_poblacion(poblacion, beta_min, beta_max, gamma_min, gamma_max, a, b, S0, I0, R0, h, indices, I_obs)
    mejores.append(max(fitnesses))
    promedios.append(sum(fitnesses) / len(fitnesses))

    #elitismo
    mejor_indice = fitnesses.index(max(fitnesses))
    mejor_individuo = poblacion[mejor_indice][:]
    
    #seccionamos usando ruleta, aplicamos cruzamiento y mutacion para generar la nueva poblacion
    seleccionados = seleccion_ruleta(poblacion, fitnesses, n_poblacion )
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
fitnesses_final = evaluar_poblacion(poblacion, beta_min, beta_max, gamma_min, gamma_max, a, b,S0, I0, R0, h, indices, I_obs )
#obtener el mejor individuo de final
mejor_indice_final = fitnesses_final.index(max(fitnesses_final))
mejor=poblacion[mejor_indice_final]

#decodificamos el mejor individuo para obtener los valores de beta y gamma
beta_final, gamma_final = decode_beta_gamma(mejor, beta_min, beta_max, gamma_min, gamma_max)
print("\nMejor solución encontrada:")
print(f"beta: {beta_final:.3f}, gamma: {gamma_final:.3f}")
print(f"R0 = {beta_final / gamma_final:.2f}")

#error relativo porcentual
error_beta = abs(beta_final - 0.3) / 0.3 * 100
error_gamma = abs(gamma_final - 0.1) / 0.1 * 100
print(f"Error relativo porcentual en beta: {error_beta:.2f}%")
print(f"Error relativo porcentual en gamma: {error_gamma:.2f}%")

#grafico del AG
graficar_modelo_AG(mejores, promedios)

#simulacion final con beta* y gamma*
f_final, jac_final = crear_modelo_sir(beta_final, gamma_final)
resultados_final = metod_taylor_segundo_orden(f_final, jac_final, a, b, [S0, I0, R0], h)
t_final = [t for t, _ in resultados_final]
S_final = [estado[0] for _, estado in resultados_final]
I_final = [estado[1] for _, estado in resultados_final]
R_final = [estado[2] for _, estado in resultados_final]

#grafico comparativo
grafico_comparacion_curvas(t_final, S_final, I_final, R_final, indices, I_obs)

'''Se busca ver como actua el algortimo AG con otros valores'''
print('-'*50)
print('Variamos los valores de L para observar los errores de cuantizacion')
L_vals = [4, 8, 16]
#para beta con p_min = 0.05 y p_max = 1
p_min = 0.05
p_max = 1

print('Los valores de errores de cuantizacion minimos son:')
for i in L_vals:
    print(f'Para un L = {i}')
    error_cuantizacion_min = (p_max - p_min)/(2**i - 1)  # Corregir: 2**i en lugar de 2*i
    print(f'Obtenemos un error = {error_cuantizacion_min:.6f} aproximadamente')
    #se observa que para un L mas grande el error de cuantizacion minima es mas precida pero hay un piso en la precision que no va a superar AG

print('-'*50)
print('Analisis de efecto de ampliar rangos de busqueda')
# Ejecutar comparación con diferentes rangos
resultados = []

# Caso 1: Rangos originales (estrechos)
caso1 = ejecutar_ag_con_rangos(0.05, 1.0, 0.01, 0.5)
resultados.append(caso1)

# Caso 2: Rangos ampliados
caso2 = ejecutar_ag_con_rangos(0.01, 2.0, 0.001, 1.0)
resultados.append(caso2)

# Caso 3: Rangos muy ampliados
caso3 = ejecutar_ag_con_rangos(0.001, 5.0, 0.0001, 2.0)
resultados.append(caso3)

for i in resultados:
    print(f"Probabilidad aproximada de encontrar solución óptima por azar: {i:.2e}, en cada caso respectivamente")
