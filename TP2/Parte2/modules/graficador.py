import matplotlib.pyplot as plt

def graficar_modelo_AG(mejores_indiv, promedios):
    plt.plot(mejores_indiv, label='Mejor fitness')
    plt.plot(promedios, label='Fitness promedio')
    plt.xlabel('Generación')
    plt.ylabel('Fitness')
    plt.legend()
    plt.title("Evolución del algoritmo genético")
    plt.show()

def grafico_comparacion_curvas(t_final, S_final, I_final, R_final, indices, I_obs):
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

def graficar_vector_I(t_values, I_values):
    plt.plot(t_values, I_values, label='I(t)')
    plt.xlabel('t')
    plt.ylabel('I')
    plt.title('Evolución de la población infectada')
    plt.legend()
    plt.show()