import matplotlib.pyplot as plt
import numpy as np

# elegi graficar cada variable (CD4, CD8 y Virus) en subplots separados en vez  de un solo gráfico,
# porque están en escalas muy diferentes, entonces en un solo gráfico el virus casi no se veía
#deje comentado el código para graficar todo en un solo gráfico por si las dudas

def graficar_resultados(simulacion, labels):
    '''Grafica los resultados de las simulaciones. Cada elemento de resultados es una lista de tuplas (t, y)'''

    tiempos = [p[0]*365 for p in simulacion] #convertimos tiempo en dias (dt=0.01, aprox 3.65 dias por paso)  -> originalmente esta en años
    CD4 = [p[1][0] for p in simulacion]
    CD8 = [p[1][1] for p in simulacion]
    V   = [p[1][2] for p in simulacion]

    fig,axs= plt.subplots(3,1,figsize=(10,8), sharex=True)

    axs[0].plot(tiempos, CD4, label=f'CD4 - {labels}', color='purple')
    axs[0].set_ylabel('CD4')
    axs[0].grid(True)

    axs[1].plot(tiempos, CD8, label=f'CD8 - {labels}', color='magenta')
    axs[1].set_ylabel('CD8')
    axs[1].grid(True)

    axs[2].plot(tiempos, V, label=f'Virus - {labels}', color='red')
    axs[2].set_xlabel('Tiempo (días)')  
    axs[2].set_ylabel('Virus')
    axs[2].grid(True)

    plt.suptitle(f"Simulación de la dinámica del VIH - {labels}") 
    plt.tight_layout()
    plt.show()
   

# Para una sola grafica:

def graf_res (simulacion, labels):
    tiempos = [p[0] for p in simulacion] 
    CD4 = [p[1][0] for p in simulacion]
    CD8 = [p[1][1] for p in simulacion]
    V   = [p[1][2]*10000 for p in simulacion] #copias/ml (multiplicamos para que se vea en misma escala que CD4 y CD8)

    plt.plot(tiempos, CD4, label=f'CD4 - {labels}', color='purple')
    plt.plot(tiempos, CD8, label=f'CD8 - {labels}', color='magenta')
    plt.plot(tiempos, V, label=f'Virus - {labels}', color='red')

    plt.xlabel('Tiempo (años)')
    plt.ylabel('Población')
    plt.title('Simulación de la dinámica del VIH')

    plt.ylim(0,1000)

    plt.legend()
    plt.grid(True)
    plt.show()
# def graficar_resultados(resultado_1, resultado_2):
#     '''Grafica los resultados de un método numérico.'''
#     fig, ax = plt.subplots()
#     tiempos = [p[0] for p in resultado_1] 
#     CD4 = [p[1][0] for p in resultado_1]
#     CD8 = [p[1][1] for p in resultado_1]
#     V   = [p[1][2]*10000 for p in resultado_1]
#     tiempos_2 = [p[0] for p in resultado_2]
#     CD4_2 = [p[1][0] for p in resultado_2]
#     CD8_2 = [p[1][1] for p in resultado_2]
#     V_2   = [p[1][2]*10000 for p in resultado_2]
#     plt.plot(tiempos, CD4, label=f'CD4 - ', color='purple')
#     plt.plot(tiempos, CD8, label=f'CD8 - ', color='magenta')
#     plt.plot(tiempos, V, label=f'Virus - ', color='red')
#     plt.plot(tiempos_2, CD4_2,  color='purple')
#     plt.plot(tiempos_2, CD8_2,  color='magenta')
#     plt.plot(tiempos_2, V_2,  color='red')
#     ax.set_xlabel('Tiempo (años)')
#     ax.set_ylabel('Población')
#     ax.set_ylim((0, 1000))
#     ax.set_xlim((0, 10))
#     ax.grid()
#     ax.legend()
#     plt.show()

