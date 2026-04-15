import random
from modules.operadores_geneticos import (evaluar_poblacion,seleccion_ruleta,cruzamiento,mutacion,decode_beta_gamma)

def correr_AG(n_poblacion, generaciones, L,beta_min, beta_max, gamma_min, gamma_max,a, b, S0, I0, R0, h, indices, I_obs):

    #población inicial aleatoria
    poblacion = [
        [random.randint(0, 1) for _ in range(L)]
        for _ in range(n_poblacion)
    ]

    mejores = []
    promedios = []

    #bucle 
    for gen in range(generaciones):
        fitnesses = evaluar_poblacion(poblacion,beta_min, beta_max, gamma_min, gamma_max,a, b, S0, I0, R0, h, indices, I_obs)
        mejores.append(max(fitnesses))
        promedios.append(sum(fitnesses)/len(fitnesses))

        #elitismo (guardar mejor individuo)
        mejor = poblacion[fitnesses.index(max(fitnesses))][:]

        #seleccionamos usando ruleta, aplicamos cruzamiento y mutacion para generar la nueva poblacion
        seleccionados = seleccion_ruleta(poblacion, fitnesses, n_poblacion)
        nueva_poblacion = []
        #cruzamiento + mutación
        for i in range(0, n_poblacion, 2):
            p1 = seleccionados[i]
            p2 = seleccionados[i+1]

            h1, h2 = cruzamiento(p1, p2)

            nueva_poblacion.append(mutacion(h1))
            nueva_poblacion.append(mutacion(h2))

        #aplicamos elitismo
        nueva_poblacion[0] = mejor

        poblacion = nueva_poblacion

    #mejor solución final
    fitness_final = evaluar_poblacion(poblacion,beta_min, beta_max, gamma_min, gamma_max,a, b, S0, I0, R0, h, indices, I_obs)

    mejor = poblacion[fitness_final.index(max(fitness_final))]
    beta_final, gamma_final = decode_beta_gamma(mejor, beta_min, beta_max, gamma_min, gamma_max)

    return beta_final, gamma_final, mejores, promedios